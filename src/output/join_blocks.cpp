/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#include <memory>
#include <thread>
#include "output.h"
#include "../util/io/temp_file.h"
#include "../data/queries.h"
#include "../output/daa/daa_write.h"
#include "output_format.h"
#include "../align/legacy/query_mapper.h"
#include "target_culling.h"
#include "../util/log_stream.h"
#include "../align/global_ranking/global_ranking.h"
#include "../legacy/util/task_queue.h"

using std::thread;
using std::unique_ptr;
using std::set;
using std::string;
using std::vector;

struct JoinFetcher
{
	static void init(const PtrVector<TempFile> &tmp_file)
	{
		for (PtrVector<TempFile>::const_iterator i = tmp_file.begin(); i != tmp_file.end(); ++i) {
			files.push_back(new InputFile(**i));
			query_ids.push_back(0);
			files.back().read(&query_ids.back(), 1);

		}
		query_last = (unsigned)-1;
	}

	static void init(const vector<string> & tmp_file_names)
	{
		for (auto file_name : tmp_file_names) {
			files.push_back(new InputFile(file_name, InputFile::NO_AUTODETECT));
			query_ids.push_back(0);
			files.back().read(&query_ids.back(), 1);
		}
		query_last = (unsigned)-1;
	}

	static void finish()
	{
		for (PtrVector<InputFile>::iterator i = files.begin(); i != files.end(); ++i)
			(*i)->close_and_delete();
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
		files[b].read(&size, 1);
		buf[b].clear();
		buf[b].resize(size);
		files[b].read(buf[b].data(), size);
		files[b].read(&query_ids[b], 1);
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
	static PtrVector<InputFile> files;
	static vector<uint32_t> query_ids;
	static unsigned query_last;
	vector<BinaryBuffer> buf;
	uint32_t query_id, unaligned_from;
};

PtrVector<InputFile> JoinFetcher::files;
vector<unsigned> JoinFetcher::query_ids;
unsigned JoinFetcher::query_last;

struct JoinWriter
{
	JoinWriter(Consumer &f):
		f_(f)
	{}
	void operator()(TextBuffer& buf)
	{
		f_.consume(buf.data(), buf.size());
		buf.clear();
	}
	Consumer &f_;
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
		std::make_heap(records.begin(), records.end(), (config.toppercent == 100.0) ? JoinRecord::cmp_evalue : JoinRecord::cmp_score);
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
		const auto pred = (config.toppercent == 100.0) ? JoinRecord::cmp_evalue : JoinRecord::cmp_score;
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
	Output::Info info = { cfg.query->seq_info(query), false, cfg.db.get(), out, {} };
	const double query_self_aln_score = flag_any(cfg.output_format->flags, Output::Flags::SELF_ALN_SCORES) ? cfg.query->self_aln_score(query) : 0.0;
	BlockJoiner joiner(buf, *cfg.db, cfg.output_format.get());
	vector<IntermediateRecord> target_hsp;
	unique_ptr<TargetCulling> culling(TargetCulling::get(cfg.max_target_seqs));

	unsigned n_target_seq = 0;
	int64_t block_idx = 0;
	OId target_oid;

	while (joiner.get(target_hsp, block_idx, target_oid, *cfg.db, cfg.output_format.get())) {
		const DictId dict_id = target_hsp.front().target_dict_id;
		const set<TaxId> rank_taxon_ids = config.taxon_k ? cfg.db->taxon_nodes().rank_taxid(cfg.db->taxids(target_oid), Rank::species) : set<TaxId>();
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
				Hsp hsp(*i, query_source_len, query_seq.index(frame).length(), tlen, cfg.output_format.get());
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

void join_worker(TaskQueue<TextBuffer, JoinWriter> *queue, const Search::Config* cfg, BitVector* ranking_db_filter_out)
{
	try {
		static std::mutex mtx;
		JoinFetcher fetcher(JoinFetcher::block_count());
		size_t n;
		TextBuffer* out;
		Statistics stat;
		const StringSet& qids = cfg->query->ids();
		//BitVector ranking_db_filter(config.global_ranking_targets > 0 ? cfg->db_seqs : 0);

		while (queue->get(n, out, fetcher) && fetcher.query_id != IntermediateRecord::FINISHED) {
			//if (!config.global_ranking_targets) stat.inc(Statistics::ALIGNED);
			stat.inc(Statistics::ALIGNED);
			size_t seek_pos;

			const char* query_name = qids[qids.check_idx(fetcher.query_id)];

			const Sequence query_seq = align_mode.query_translated ? cfg->query->source_seqs()[fetcher.query_id] : cfg->query->seqs()[fetcher.query_id];

			if (*cfg->output_format != OutputFormat::daa && config.report_unaligned != 0) {
				for (unsigned i = fetcher.unaligned_from; i < fetcher.query_id; ++i) {
					Output::Info info{ cfg->query->seq_info(i), true, cfg->db.get(), *out, {} };
					cfg->output_format->print_query_intro(info);
					cfg->output_format->print_query_epilog(info);
				}
			}

			unique_ptr<OutputFormat> f(cfg->output_format->clone());

			Output::Info info{ cfg->query->seq_info(fetcher.query_id), false, cfg->db.get(), *out, {} };
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
	catch (std::exception& e) {
		exit_with_error(e);
	}
}

void join_blocks(int64_t ref_blocks, Consumer &master_out, const PtrVector<TempFile> &tmp_file, Search::Config& cfg, SequenceFile &db_file,
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

	unique_ptr<TempFile> merged_query_list;
	/*if (config.global_ranking_targets)
		merged_query_list.reset(new TempFile());
	JoinWriter writer(config.global_ranking_targets ? *merged_query_list : master_out);*/
	JoinWriter writer(master_out);
	TaskQueue<TextBuffer, JoinWriter> queue(3 * config.threads_, writer);
	vector<thread> threads;
	//BitVector ranking_db_filter(config.global_ranking_targets > 0 ? cfg.db_seqs : 0);
	for (int i = 0; i < config.threads_; ++i)
		threads.emplace_back(join_worker, &queue, &cfg, nullptr); // &ranking_db_filter);
	for (auto &t : threads)
		t.join();
	JoinFetcher::finish();
	if (*cfg.output_format != OutputFormat::daa && config.report_unaligned != 0) {
		TextBuffer out;
		for (BlockId i = JoinFetcher::query_last + 1; i < query_ids.size(); ++i) {
			Output::Info info{ cfg.query->seq_info(i), true, cfg.db.get(), out, {} };
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