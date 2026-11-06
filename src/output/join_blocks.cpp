/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/
// SPDX-License-Identifier: GPL-3.0-or-later

#include <memory>
#include <thread>
#include "output.h"
#include "legacy/daa/daa_write.h"
#include "output_format.h"
#include "target_culling.h"
#include "util/log_stream.h"
#include "legacy/util/task_queue.h"
#include "data/taxonomy_nodes.h"
#include "util/parallel/simple_thread_pool.h"

using std::thread;
using std::unique_ptr;
using std::set;
using std::string;
using std::vector;

struct JoinFetcher
{
	static void init(const vector<File*> &tmp_file)
	{
		for (vector<File*>::const_iterator i = tmp_file.begin(); i != tmp_file.end(); ++i) {
			(*i)->rewind();
			files.push_back(*i);
			query_ids.push_back(0);
			files.back()->read(&query_ids.back(), sizeof(uint32_t));

		}
		query_last = (unsigned)-1;
	}

	static void init(const vector<string> & tmp_file_names)
	{
		for (auto file_name : tmp_file_names) {
			files.push_back(new File(file_name, "rb"));
			query_ids.push_back(0);
			files.back()->read(&query_ids.back(), sizeof(uint32_t));
		}
		query_last = (unsigned)-1;
	}

	static void finish()
	{
		for (vector<File*>::iterator i = files.begin(); i != files.end(); ++i) {
			(*i)->close();
			delete* i;
		}
		files.clear();
		query_ids.clear();
	}
	static uint32_t next()
	{
		return *std::min_element(query_ids.begin(), query_ids.end());
	}
	static size_t block_count() {
		return files.size();
	}
	void fetch(unsigned b)
	{
		uint32_t size;
		files[b]->read(&size, sizeof(uint32_t));
		buf[b].clear();
		buf[b].resize(size);
		files[b]->read(buf[b].data(), size);
		files[b]->read(&query_ids[b], sizeof(uint32_t));
	}
	JoinFetcher(size_t blocks):
		buf(blocks)
	{}
	bool operator()()
	{
		query_id = next();
		unaligned_from = query_last + 1;
		query_last = query_id;
		for (unsigned i = 0; i < buf.size(); ++i)
			if (query_ids[i] == query_id && query_id != IntermediateRecord::FINISHED)
				fetch(i);
			else
				buf[i].clear();
		return next() != IntermediateRecord::FINISHED;
	}
	static vector<File*> files;
	static vector<uint32_t> query_ids;
	static unsigned query_last;
	vector<BinaryBuffer> buf;
	uint32_t query_id, unaligned_from;
};

vector<File*> JoinFetcher::files;
vector<unsigned> JoinFetcher::query_ids;
unsigned JoinFetcher::query_last;

struct JoinWriter
{
	JoinWriter(File &f):
		f_(f)
	{}
	void operator()(TextBuffer& buf)
	{
		f_.write(buf.data(), buf.size());
		buf.clear();
	}
	File &f_;
};

struct JoinRecord
{

	static bool cmp_evalue(const JoinRecord& lhs, const JoinRecord &rhs)
	{
		return rhs.same_subject_ || (!rhs.same_subject_ && (lhs.info_.evalue > rhs.info_.evalue || (lhs.info_.evalue == rhs.info_.evalue && cmp_score(lhs, rhs))));
	}

	static bool cmp_score(const JoinRecord& lhs, const JoinRecord& rhs)
	{
		return rhs.same_subject_ || (!rhs.same_subject_ && (lhs.info_.score < rhs.info_.score || (lhs.info_.score == rhs.info_.score && rhs.db_precedence(lhs))));
	}

	bool db_precedence(const JoinRecord &rhs) const
	{
		return info_.target_oid < rhs.info_.target_oid;
	}

	JoinRecord(int64_t ref_block, DictId subject, BinaryBuffer::Iterator &it, const SequenceFile& db, const OutputFormat* output_format):
		block_(ref_block)
	{
		info_.read(it, output_format);
		same_subject_ = info_.target_dict_id == subject;
		if (*output_format != OutputFormat::daa)
			info_.target_oid = db.oid(info_.target_dict_id, ref_block);
	}

	static bool push_next(int64_t block, DictId subject, BinaryBuffer::Iterator &it, vector<JoinRecord> &v, const SequenceFile& db, const OutputFormat* output_format)
	{
		if (it.good()) {
			v.push_back(JoinRecord(block, subject, it, db, output_format));
			return true;
		}
		else
			return false;
	}

	int64_t block_;
	bool same_subject_;
	IntermediateRecord info_;

};

struct BlockJoiner
{
	BlockJoiner(vector<BinaryBuffer> &buf, const SequenceFile& db, const OutputFormat* output_format)
	{
		for (auto i = buf.begin(); i != buf.end(); ++i) {
			it.push_back(i->begin());
			JoinRecord::push_next(i - buf.begin(), std::numeric_limits<unsigned>::max(), it.back(), records, db, output_format);
		}
		//std::make_heap(records.begin(), records.end(), (config.toppercent == 100.0 && config.global_ranking_targets == 0) ? JoinRecord::cmp_evalue : JoinRecord::cmp_score);
		std::make_heap(records.begin(), records.end(), config.toppercent.blank() ? JoinRecord::cmp_evalue : JoinRecord::cmp_score);
	}
	bool get(vector<IntermediateRecord> &target_hsp, int64_t& block_idx, OId& target_oid, const SequenceFile& db, const OutputFormat* output_format)
	{
		if (records.empty())
			return false;
		const JoinRecord &first = records.front();
		const int64_t block = first.block_;
		block_idx = block;
		target_oid = first.info_.target_oid;
		const DictId subject = first.info_.target_dict_id;
		target_hsp.clear();
		//const auto pred = (config.toppercent == 100.0 && config.global_ranking_targets == 0) ? JoinRecord::cmp_evalue : JoinRecord::cmp_score;
		const auto pred = config.toppercent.blank() ? JoinRecord::cmp_evalue : JoinRecord::cmp_score;
		do {
			const JoinRecord &next = records.front();
			if (next.block_ != block || next.info_.target_dict_id != subject)
				return true;
			target_hsp.push_back(next.info_);
			std::pop_heap(records.begin(), records.end(), pred);
			records.pop_back();
			if (JoinRecord::push_next(block, subject, it[block], records, db, output_format))
				std::push_heap(records.begin(), records.end(), pred);
		} while (!records.empty());
		return true;
	}
	vector<JoinRecord> records;
	vector<BinaryBuffer::Iterator> it;
};

void join_query(
	vector<BinaryBuffer> &buf,
	TextBuffer &out,
	Statistics &statistics,
	unsigned query,
	const char *query_name,
	unsigned query_source_len,
	OutputFormat &f,
	const Search::Config &cfg)
{
	TranslatedSequence query_seq(cfg.query->translated(query));
	Output::Info info = { cfg.query->seq_info(query), false, cfg.db.get(), out, Util::Seq::AccessionParsing(), cfg.db->sequence_count(), cfg.db->letters() };
	const double query_self_aln_score = flag_any(cfg.output_format->flags, Output::Flags::SELF_ALN_SCORES) ? cfg.query->self_aln_score(query) : 0.0;
	BlockJoiner joiner(buf, *cfg.db, cfg.output_format.get());
	vector<IntermediateRecord> target_hsp;
	unique_ptr<TargetCulling> culling(TargetCulling::get(cfg.max_target_seqs));

	unsigned n_target_seq = 0;
	int64_t block_idx = 0;
	OId target_oid;

	while (joiner.get(target_hsp, block_idx, target_oid, *cfg.db, cfg.output_format.get())) {
		const DictId dict_id = target_hsp.front().target_dict_id;
		const set<TaxId> rank_taxon_ids = config.taxon_k ? cfg.db->rank_taxid(cfg.db->taxids(target_oid), Rank::species) : set<TaxId>();
		const int c = culling->cull(target_hsp, rank_taxon_ids);
		const double target_self_aln_score = flag_any(cfg.output_format->flags, Output::Flags::SELF_ALN_SCORES) ? cfg.db->dict_self_aln_score(dict_id, block_idx) : 0.0;
		if (c == TargetCulling::FINISHED)
			break;
		else if (c == TargetCulling::NEXT)
			continue;

		unsigned hsp_num = 0;
		for (vector<IntermediateRecord>::const_iterator i = target_hsp.begin(); i != target_hsp.end(); ++i, ++hsp_num) {
			if (f == OutputFormat::daa)
				write_daa_record(out, *i);
			/*else if (config.global_ranking_targets > 0)
				Extension::GlobalRanking::write_merged_query_list(*i, out, ranking_db_filter, statistics);*/
			else {
				const Loc tlen = cfg.db->dict_len(dict_id, block_idx);
				const unsigned frame = i->frame(query_source_len, align_mode.mode);

#ifdef WITH_DNA
                Hsp hsp = Hsp(*i, query_source_len, query_seq.index(frame).length(), tlen, cfg.output_format.get(), cfg.score_builder.get());
#else
				Hsp hsp = Hsp(*i, query_source_len, query_seq.index(frame).length(), tlen, cfg.output_format.get());
#endif

				f.print_match(HspContext(hsp,
					query,
					cfg.query->block_id2oid(query),
					query_seq,
					query_name,
					target_oid,
					tlen,
					cfg.db->dict_title(dict_id, block_idx).c_str(),
					n_target_seq,
					hsp_num,
					flag_any(f.flags, Output::Flags::TARGET_SEQS) ? Sequence(cfg.db->dict_seq(dict_id, block_idx)) : Sequence(),
					0,
					query_self_aln_score,
					target_self_aln_score).parse(cfg.output_format.get()), info);
			}
		}

		culling->add(target_hsp, rank_taxon_ids);
		++n_target_seq;
		//if (!config.global_ranking_targets) {
			statistics.inc(Statistics::PAIRWISE);
			statistics.inc(Statistics::MATCHES, hsp_num);
		//}
	}
}

void join_worker(const std::atomic<bool>& stop, TaskQueue<TextBuffer, JoinWriter>* queue, const Search::Config* cfg, BitVector* ranking_db_filter_out)
{
	static std::mutex mtx;
	JoinFetcher fetcher(JoinFetcher::block_count());
	size_t n;
	TextBuffer* out;
	Statistics stat;
	const StringSet& qids = cfg->query->ids();
	//BitVector ranking_db_filter(config.global_ranking_targets > 0 ? cfg->db_seqs : 0);

	while (!stop.load(std::memory_order_relaxed) && queue->get(n, out, fetcher) && fetcher.query_id != IntermediateRecord::FINISHED) {
		//if (!config.global_ranking_targets) stat.inc(Statistics::ALIGNED);
		stat.inc(Statistics::ALIGNED);
		size_t seek_pos;

		const char* query_name = qids[qids.check_idx(fetcher.query_id)];

		const Sequence query_seq = align_mode.query_translated ? cfg->query->source_seqs()[fetcher.query_id] : cfg->query->seqs()[fetcher.query_id];

		if (n > 0 && cfg->output_format->query_separator != '\0')
			out->write(cfg->output_format->query_separator);
		if (*cfg->output_format != OutputFormat::daa && config.report_unaligned != 0) {
			for (unsigned i = fetcher.unaligned_from; i < fetcher.query_id; ++i) {
				Output::Info info{ cfg->query->seq_info(i), true, cfg->db.get(), *out, Util::Seq::AccessionParsing(), cfg->db->sequence_count(), cfg->db->letters() };
				cfg->output_format->print_query_intro(info);
				cfg->output_format->print_query_epilog(info);
			}
		}

		unique_ptr<OutputFormat> f(cfg->output_format->clone());

		Output::Info info{ cfg->query->seq_info(fetcher.query_id), false, cfg->db.get(), *out, Util::Seq::AccessionParsing(), cfg->db->sequence_count(), cfg->db->letters() };
		if (*f == OutputFormat::daa)
			seek_pos = write_daa_query_record(*out, query_name, query_seq);
		/*else if (config.global_ranking_targets)
			seek_pos = Extension::GlobalRanking::write_merged_query_list_intro(fetcher.query_id, *out);*/
		else
			f->print_query_intro(info);

		join_query(fetcher.buf, *out, stat, fetcher.query_id, query_name, (unsigned)query_seq.length(), *f, *cfg); // ranking_db_filter);

		if (*f == OutputFormat::daa)
			finish_daa_query_record(*out, seek_pos);
		/*else if (config.global_ranking_targets)
			Extension::GlobalRanking::finish_merged_query_list(*out, seek_pos);*/
		else
			f->print_query_epilog(info);
		queue->push(n);
	}

	statistics += stat;
	/*if (config.global_ranking_targets) {
		std::lock_guard<std::mutex> lock(mtx);
		(*ranking_db_filter_out) |= ranking_db_filter;
	}*/
}

struct JoinBatch {
	static constexpr size_t MAX_Q = 64, MAX_BYTES = 256 << 10;
	struct Slice { size_t off; uint32_t len; };
	BinaryBuffer raw; vector<Slice> slices;
	uint32_t qids[MAX_Q], unaligned_from[MAX_Q], count = 0; size_t blocks;
	JoinBatch(size_t b) : slices(MAX_Q * b), blocks(b) {}
	bool operator()() {
		raw.clear(); count = 0;
		while (count < MAX_Q && (count == 0 || raw.size() < MAX_BYTES)) {
			const uint32_t q = JoinFetcher::next(); if (q == IntermediateRecord::FINISHED) break;
			qids[count] = q; unaligned_from[count] = JoinFetcher::query_last + 1; JoinFetcher::query_last = q;
			Slice* row = &slices[count * blocks];
			for (size_t b = 0; b < blocks; ++b) if (JoinFetcher::query_ids[b] == q) {
				uint32_t sz; JoinFetcher::files[b]->read(&sz, sizeof(uint32_t));
				const size_t off = raw.size(); raw.resize(off + sz); JoinFetcher::files[b]->read(raw.data() + off, sz);
				JoinFetcher::files[b]->read(&JoinFetcher::query_ids[b], sizeof(uint32_t)); row[b] = { off, sz };
			} else row[b] = { 0, 0 };
			++count;
		}
		return JoinFetcher::next() != IntermediateRecord::FINISHED;
	}
};
static void join_worker_fast(const std::atomic<bool>& stop, TaskQueue<TextBuffer, JoinWriter>* queue, const Search::Config* cfg) {
	const size_t blocks = JoinFetcher::block_count(); JoinBatch batch(blocks);
	vector<JoinRecord> records; vector<BinaryBuffer::Iterator> it; vector<IntermediateRecord> target_hsp; vector<Letter> target_seq;
	HspContext ctx; unique_ptr<OutputFormat> fmt(cfg->output_format->clone()); OutputFormat& f = *fmt;
	const bool rep_unal = f != OutputFormat::daa && cfg->output_format->report_unaligned();
	const auto pred = config.toppercent.blank() ? JoinRecord::cmp_evalue : JoinRecord::cmp_score;
	const auto db_seqs = cfg->db->sequence_count(), db_letters = cfg->db->letters();
	size_t n; TextBuffer* out; Statistics stat; const StringSet& qids = cfg->query->ids();
	while (!stop.load(std::memory_order_relaxed) && queue->get(n, out, batch) && batch.count > 0) {
		stat.inc(Statistics::ALIGNED, batch.count);
		for (uint32_t k = 0; k < batch.count; ++k) {
			const uint32_t qid = batch.qids[k]; const char* query_name = qids[qids.check_idx(qid)];
			const Sequence qsrc = align_mode.query_translated ? cfg->query->source_seqs()[qid] : cfg->query->seqs()[qid];
			const unsigned qlen = (unsigned)qsrc.length();
			if ((n > 0 || k > 0) && cfg->output_format->query_separator != '\0') out->write(cfg->output_format->query_separator);
			if (rep_unal) for (unsigned i = batch.unaligned_from[k]; i < qid; ++i) {
				Output::Info ui{ cfg->query->seq_info(i), true, cfg->db.get(), *out, Util::Seq::AccessionParsing(), db_seqs, db_letters };
				cfg->output_format->print_query_intro(ui); cfg->output_format->print_query_epilog(ui);
			}
			if (f == OutputFormat::taxon) { static_cast<TaxonFormat*>(&f)->taxid = 0; static_cast<TaxonFormat*>(&f)->evalue = std::numeric_limits<double>::max(); }
			Output::Info info{ cfg->query->seq_info(qid), false, cfg->db.get(), *out, Util::Seq::AccessionParsing(), db_seqs, db_letters };
			size_t seek_pos = f == OutputFormat::daa ? write_daa_query_record(*out, query_name, qsrc) : (f.print_query_intro(info), (size_t)0);
			TranslatedSequence query_seq(cfg->query->translated(qid)); const JoinBatch::Slice* slices = &batch.slices[k * blocks]; records.clear(); it.clear();
			for (size_t b = 0; b < blocks; ++b) {
				it.emplace_back(batch.raw.cbegin() + slices[b].off, batch.raw.cbegin() + slices[b].off + slices[b].len);
				JoinRecord::push_next((int64_t)b, std::numeric_limits<unsigned>::max(), it.back(), records, *cfg->db, cfg->output_format.get());
			}
			std::make_heap(records.begin(), records.end(), pred);
			GlobalCulling gc(cfg->max_target_seqs);
			unique_ptr<TargetCulling> rc(config.query_range_culling ? TargetCulling::get(cfg->max_target_seqs) : nullptr);
			TargetCulling* culling = rc ? rc.get() : &gc;
			ctx.query = query_seq; ctx.query_title.assign(query_name); ctx.query_id = qid; ctx.query_oid = cfg->query->block_id2oid(qid);
			ctx.query_len = query_seq.source().length();
			ctx.query_self_aln_score = flag_any(cfg->output_format->flags, Output::Flags::SELF_ALN_SCORES) ? cfg->query->self_aln_score(qid) : 0.0;
			unsigned n_target_seq = 0;
			while (!records.empty()) {
				const int64_t bidx = records.front().block_; const OId toid = records.front().info_.target_oid; const DictId did = records.front().info_.target_dict_id;
				target_hsp.clear();
				do {
					if (records.front().block_ != bidx || records.front().info_.target_dict_id != did) break;
					std::pop_heap(records.begin(), records.end(), pred);
					target_hsp.push_back(std::move(records.back().info_)); records.pop_back();
					if (JoinRecord::push_next(bidx, did, it[bidx], records, *cfg->db, cfg->output_format.get())) std::push_heap(records.begin(), records.end(), pred);
				} while (!records.empty());
				const set<TaxId> rtax = config.taxon_k ? cfg->db->rank_taxid(cfg->db->taxids(toid), Rank::species) : set<TaxId>();
				const int c = culling->cull(target_hsp, rtax);
				if (c == TargetCulling::FINISHED) break; else if (c == TargetCulling::NEXT) continue;
				if (f != OutputFormat::daa) {
					ctx.subject_oid = toid; ctx.subject_len = cfg->db->dict_len(did, bidx); ctx.hit_num = n_target_seq;
					ctx.target_title.assign(cfg->db->dict_title_c_str(did, bidx));
					ctx.target_self_aln_score = flag_any(cfg->output_format->flags, Output::Flags::SELF_ALN_SCORES) ? cfg->db->dict_self_aln_score(did, bidx) : 0.0;
					ctx.subject_seq = flag_any(f.flags, Output::Flags::TARGET_SEQS) ? (target_seq = cfg->db->dict_seq(did, bidx), Sequence(target_seq)) : Sequence();
				}
				unsigned hsp_num = 0;
				for (auto i = target_hsp.cbegin(); i != target_hsp.cend(); ++i, ++hsp_num)
					if (f == OutputFormat::daa) write_daa_record(*out, *i);
					else {
						const unsigned frame = i->frame(qlen, align_mode.mode);
#ifdef WITH_DNA
						ctx.set_hsp(Hsp(*i, qlen, query_seq.index(frame).length(), ctx.subject_len, cfg->output_format.get(), cfg->score_builder.get()));
#else
						ctx.set_hsp(Hsp(*i, qlen, query_seq.index(frame).length(), ctx.subject_len, cfg->output_format.get()));
#endif
						ctx.hsp_num = hsp_num; f.print_match(ctx.parse(cfg->output_format.get()), info);
					}
				culling->add(target_hsp, rtax);
				++n_target_seq; stat.inc(Statistics::PAIRWISE); stat.inc(Statistics::MATCHES, hsp_num);
			}
			if (f == OutputFormat::daa) finish_daa_query_record(*out, seek_pos); else f.print_query_epilog(info);
		}
		queue->push(n);
	}
	statistics += stat;
}
void join_blocks(int64_t ref_blocks, File &master_out, const vector<File*> &tmp_file, Search::Config& cfg, SequenceFile &db_file,
	const vector<string> tmp_file_names)
{
	if (*cfg.output_format != OutputFormat::daa)
		cfg.db->init_random_access(cfg.current_query_block, config.multiprocessing ? tmp_file_names.size() : tmp_file.size());
	TaskTimer timer("Joining output blocks");

	if (tmp_file_names.size() > 0) {
		JoinFetcher::init(tmp_file_names);
	} else {
		JoinFetcher::init(tmp_file);
	}

	const StringSet& query_ids = cfg.query->ids();

	//unique_ptr<TempFile> merged_query_list;
	/*if (config.global_ranking_targets)
		merged_query_list.reset(new TempFile());
	JoinWriter writer(config.global_ranking_targets ? *merged_query_list : master_out);*/
	JoinWriter writer(master_out);
	TaskQueue<TextBuffer, JoinWriter> queue(3 * config.threads_, writer);
	SimpleThreadPool pool;
	//BitVector ranking_db_filter(config.global_ranking_targets > 0 ? cfg.db_seqs : 0);
	for (int i = 0; i < config.threads_; ++i)
		pool.spawn(join_worker_fast, &queue, &cfg);
	pool.join_all();
	JoinFetcher::finish();
	if (*cfg.output_format != OutputFormat::daa && cfg.output_format->report_unaligned()) {
		TextBuffer out;
		for (BlockId i = JoinFetcher::query_last + 1; i < query_ids.size(); ++i) {
			Output::Info info{ cfg.query->seq_info(i), true, cfg.db.get(), out, Util::Seq::AccessionParsing(), cfg.db->sequence_count(), cfg.db->letters() };
			cfg.output_format->print_query_intro(info);
			cfg.output_format->print_query_epilog(info);
		}
		writer(out);
	}

	/*if (config.global_ranking_targets)
		Extension::GlobalRanking::extend(db_file, *merged_query_list, ranking_db_filter, cfg, master_out);*/
	if (*cfg.output_format != OutputFormat::daa)
		cfg.db->end_random_access();
}