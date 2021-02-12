/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
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
#include "output.h"
#include "../util/io/temp_file.h"
#include "../data/queries.h"
#include "../output/daa_write.h"
#include "output_format.h"
#include "../align/legacy/query_mapper.h"
#include "target_culling.h"
#include "../data/ref_dictionary.h"
#include "../util/log_stream.h"
#include "../align/global_ranking/global_ranking.h"

using namespace std;

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
	void fetch(unsigned b)
	{
		uint32_t size;
		files[b].read(&size, 1);
		buf[b].clear();
		buf[b].resize(size);
		files[b].read(buf[b].data(), size);
		files[b].read(&query_ids[b], 1);
	}
	JoinFetcher():
		buf(current_ref_block)
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
		f_.consume(buf.get_begin(), buf.size());
		buf.clear();
	}
	Consumer &f_;
};

struct Join_record
{

	static bool cmp_evalue(const Join_record& lhs, const Join_record &rhs)
	{
		return rhs.same_subject_ || (!rhs.same_subject_ && (lhs.info_.evalue > rhs.info_.evalue || (lhs.info_.evalue == rhs.info_.evalue && cmp_score(lhs, rhs))));
	}

	static bool cmp_score(const Join_record& lhs, const Join_record& rhs)
	{
		return rhs.same_subject_ || (!rhs.same_subject_ && (lhs.info_.score < rhs.info_.score || (lhs.info_.score == rhs.info_.score && rhs.db_precedence(lhs))));
	}

	bool db_precedence(const Join_record &rhs) const
	{
		return block_ < rhs.block_ || (block_ == rhs.block_ && info_.subject_dict_id < rhs.info_.subject_dict_id);
	}

	Join_record(unsigned ref_block, unsigned subject, BinaryBuffer::Iterator &it):
		block_(ref_block)
	{
		info_.read(it);
		same_subject_ = info_.subject_dict_id == subject;
	}

	static bool push_next(unsigned block, unsigned subject, BinaryBuffer::Iterator &it, vector<Join_record> &v)
	{
		if (it.good()) {
			v.push_back(Join_record(block, subject, it));
			return true;
		}
		else
			return false;
	}

	unsigned block_;
	bool same_subject_;
	IntermediateRecord info_;

};

struct BlockJoiner
{
	BlockJoiner(vector<BinaryBuffer> &buf)
	{
		for (unsigned i = 0; i < current_ref_block; ++i) {
			it.push_back(buf[i].begin());
			Join_record::push_next(i, std::numeric_limits<unsigned>::max(), it.back(), records);
		}
		std::make_heap(records.begin(), records.end(), (config.toppercent == 100.0 && config.global_ranking_targets == 0) ? Join_record::cmp_evalue : Join_record::cmp_score);
	}
	bool get(vector<IntermediateRecord> &target_hsp, unsigned & block_idx)
	{
		if (records.empty())
			return false;
		const Join_record &first = records.front();
		const unsigned block = first.block_;
		block_idx = block;
		const unsigned subject = first.info_.subject_dict_id;
		target_hsp.clear();
		const auto pred = (config.toppercent == 100.0 && config.global_ranking_targets == 0) ? Join_record::cmp_evalue : Join_record::cmp_score;
		do {
			const Join_record &next = records.front();
			if (next.block_ != block || next.info_.subject_dict_id != subject)
				return true;
			target_hsp.push_back(next.info_);
			std::pop_heap(records.begin(), records.end(), pred);
			records.pop_back();
			if (Join_record::push_next(block, subject, it[block], records))
				std::push_heap(records.begin(), records.end(), pred);
		} while (!records.empty());
		return true;
	}
	vector<Join_record> records;
	vector<BinaryBuffer::Iterator> it;
};

void join_query(
	vector<BinaryBuffer> &buf,
	TextBuffer &out,
	Statistics &statistics,
	unsigned query,
	const char *query_name,
	unsigned query_source_len,
	Output_format &f,
	const Metadata &metadata,
	BitVector& ranking_db_filter)
{
	ReferenceDictionary& dict = ReferenceDictionary::get();
	TranslatedSequence query_seq(get_translated_query(query));
	BlockJoiner joiner(buf);
	vector<IntermediateRecord> target_hsp;
	unique_ptr<TargetCulling> culling(TargetCulling::get());

	unsigned n_target_seq = 0;
	unsigned block_idx = 0;

	while (joiner.get(target_hsp, block_idx)) {
		ReferenceDictionary * dict_ptr;
		if (config.multiprocessing) {
			dict_ptr = & ReferenceDictionary::get(block_idx);
		} else {
			dict_ptr = & dict;
		}

		const set<unsigned> rank_taxon_ids = config.taxon_k ? metadata.taxon_nodes->rank_taxid((*metadata.taxon_list)[dict_ptr->database_id(target_hsp.front().subject_dict_id)], Rank::species) : set<unsigned>();
		const int c = culling->cull(target_hsp, rank_taxon_ids);
		if (c == TargetCulling::FINISHED)
			break;
		else if (c == TargetCulling::NEXT)
			continue;

		unsigned hsp_num = 0;
		for (vector<IntermediateRecord>::const_iterator i = target_hsp.begin(); i != target_hsp.end(); ++i, ++hsp_num) {
			if (f == Output_format::daa)
				write_daa_record(out, *i);
			else if (config.global_ranking_targets > 0)
				Extension::GlobalRanking::write_merged_query_list(*i, dict, out, ranking_db_filter, statistics);
			else {
				Hsp hsp(*i, query_source_len);
				f.print_match(Hsp_context(hsp,
					query,
					query_seq,
					query_name,
					dict_ptr->check_id(i->subject_dict_id),
					dict_ptr->database_id(i->subject_dict_id),
					dict_ptr->name(i->subject_dict_id),
					dict_ptr->length(i->subject_dict_id),
					n_target_seq,
					hsp_num,
					config.use_lazy_dict ? dict_ptr->seq(i->subject_dict_id) : Sequence()
					).parse(), metadata, out);
			}
		}

		culling->add(target_hsp, rank_taxon_ids);
		++n_target_seq;
		if (!config.global_ranking_targets) {
			statistics.inc(Statistics::PAIRWISE);
			statistics.inc(Statistics::MATCHES, hsp_num);
		}
	}
}

void join_worker(Task_queue<TextBuffer, JoinWriter> *queue, const Parameters *params, const Metadata *metadata, BitVector* ranking_db_filter_out)
{
	static std::mutex mtx;
	JoinFetcher fetcher;
	size_t n;
	TextBuffer *out;
	Statistics stat;
	const String_set<char, 0>& qids = query_ids::get();
	BitVector ranking_db_filter(config.global_ranking_targets > 0 ? params->db_seqs : 0);

	while (queue->get(n, out, fetcher) && fetcher.query_id != IntermediateRecord::FINISHED) {
		if(!config.global_ranking_targets) stat.inc(Statistics::ALIGNED);
		size_t seek_pos;

		const char * query_name = qids[qids.check_idx(fetcher.query_id)];

		const Sequence query_seq = align_mode.query_translated ? query_source_seqs::get()[fetcher.query_id] : query_seqs::get()[fetcher.query_id];

		if (*output_format != Output_format::daa && config.report_unaligned != 0) {
			for (unsigned i = fetcher.unaligned_from; i < fetcher.query_id; ++i) {
				output_format->print_query_intro(i, query_ids::get()[i], get_source_query_len(i), *out, true);
				output_format->print_query_epilog(*out, query_ids::get()[i], true, *params);
			}
		}

		unique_ptr<Output_format> f(output_format->clone());

		if (*f == Output_format::daa)
			seek_pos = write_daa_query_record(*out, query_name, query_seq);
		else if (config.global_ranking_targets)
			seek_pos = Extension::GlobalRanking::write_merged_query_list_intro(fetcher.query_id, *out);
		else
			f->print_query_intro(fetcher.query_id, query_name, (unsigned)query_seq.length(), *out, false);

		join_query(fetcher.buf, *out, stat, fetcher.query_id, query_name, (unsigned)query_seq.length(), *f, *metadata, ranking_db_filter);

		if (*f == Output_format::daa)
			finish_daa_query_record(*out, seek_pos);
		else if (config.global_ranking_targets)
			Extension::GlobalRanking::finish_merged_query_list(*out, seek_pos);
		else
			f->print_query_epilog(*out, query_name, false, *params);
		queue->push(n);
	}

	statistics += stat;
	if (config.global_ranking_targets) {
		std::lock_guard<std::mutex> lock(mtx);
		(*ranking_db_filter_out) |= ranking_db_filter;
	}
}

void join_blocks(unsigned ref_blocks, Consumer &master_out, const PtrVector<TempFile> &tmp_file, const Parameters &params, const Metadata &metadata, DatabaseFile &db_file,
	const vector<string> tmp_file_names)
{
	//ReferenceDictionary::get().init_rev_map();
	task_timer timer("Building reference dictionary", 3);
	if (config.use_lazy_dict)
		ReferenceDictionary::get().build_lazy_dict(db_file);
	timer.go("Joining output blocks");

	if (tmp_file_names.size() > 0) {
		JoinFetcher::init(tmp_file_names);
	} else {
		JoinFetcher::init(tmp_file);
	}

	unique_ptr<TempFile> merged_query_list;
	if (config.global_ranking_targets)
		merged_query_list.reset(new TempFile());
	JoinWriter writer(config.global_ranking_targets ? *merged_query_list : master_out);
	Task_queue<TextBuffer, JoinWriter> queue(3 * config.threads_, writer);
	vector<thread> threads;
	BitVector ranking_db_filter(config.global_ranking_targets > 0 ? params.db_seqs : 0);
	for (unsigned i = 0; i < config.threads_; ++i)
		threads.emplace_back(join_worker, &queue, &params, &metadata, &ranking_db_filter);
	for (auto &t : threads)
		t.join();
	JoinFetcher::finish();
	if (*output_format != Output_format::daa && config.report_unaligned != 0) {
		TextBuffer out;
		for (unsigned i = JoinFetcher::query_last + 1; i < query_ids::get().get_length(); ++i) {
			output_format->print_query_intro(i, query_ids::get()[i], get_source_query_len(i), out, true);
			output_format->print_query_epilog(out, query_ids::get()[i], true, params);
		}
		writer(out);
	}

	if (config.use_lazy_dict) {
		timer.go("Deallocating dictionary");
		delete ref_seqs::data_;
		ref_seqs::data_ = NULL;
	}

	if (config.global_ranking_targets)
		Extension::GlobalRanking::extend(db_file, *merged_query_list, ranking_db_filter, params, metadata, master_out);
}