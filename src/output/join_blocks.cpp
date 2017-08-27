/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include "output.h"
#include "../util/temp_file.h"
#include "../data/queries.h"
#include "../output/daa_write.h"
#include "output_format.h"
#include "../align/query_mapper.h"

struct Join_fetcher
{
	static void init(const vector<Temp_file> &tmp_file)
	{
		for (vector<Temp_file>::const_iterator i = tmp_file.begin(); i != tmp_file.end(); ++i) {
			files.push_back(Input_stream(*i));
			query_ids.push_back(0);
			files.back().read(&query_ids.back(), 1);
		}
		query_last = (unsigned)-1;
	}
	static void finish()
	{
		for (vector<Input_stream>::iterator i = files.begin(); i != files.end(); ++i)
			i->close_and_delete();
		files.clear();
		query_ids.clear();
	}
	static unsigned next()
	{
		return *std::min_element(query_ids.begin(), query_ids.end());
	}
	void fetch(unsigned b)
	{
		unsigned size;
		files[b].read(&size, 1);
		buf[b].clear();
		buf[b].resize(size);
		files[b].read(buf[b].data(), size);
		files[b].read(&query_ids[b], 1);
	}
	Join_fetcher():
		buf(current_ref_block)
	{}
	bool operator()()
	{
		query_id = next();
		unaligned_from = query_last + 1;
		query_last = query_id;
		for (unsigned i = 0; i < buf.size(); ++i)
			if (query_ids[i] == query_id && query_id != Intermediate_record::finished)
				fetch(i);
			else
				buf[i].clear();
		return next() != Intermediate_record::finished;
	}
	static vector<Input_stream> files;
	static vector<unsigned> query_ids;
	static unsigned query_last;
	vector<Binary_buffer> buf;
	unsigned query_id, unaligned_from;
};

vector<Input_stream> Join_fetcher::files;
vector<unsigned> Join_fetcher::query_ids;
unsigned Join_fetcher::query_last;

struct Join_writer
{
	Join_writer(Output_stream &f):
		f_(f)
	{}
	void operator()(Text_buffer& buf)
	{
		f_.write(buf.get_begin(), buf.size());
		buf.clear();
	}
	Output_stream &f_;
};

struct Join_record
{

	bool operator<(const Join_record &rhs) const
	{
		return rhs.same_subject_ ||  (!rhs.same_subject_ && (info_.score < rhs.info_.score || (info_.score == rhs.info_.score && rhs.db_precedence(*this))));
	}

	bool db_precedence(const Join_record &rhs) const
	{
		return block_ < rhs.block_ || (block_ == rhs.block_ && info_.subject_id < rhs.info_.subject_id);
	}

	Join_record(unsigned ref_block, unsigned subject, Binary_buffer::Iterator &it):
		block_(ref_block)
	{
		info_.read(it);
		same_subject_ = info_.subject_id == subject;
	}

	static bool push_next(unsigned block, unsigned subject, Binary_buffer::Iterator &it, vector<Join_record> &v)
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
	Intermediate_record info_;

};

struct Join_records
{
	Join_records(vector<Binary_buffer> &buf)
	{
		for (unsigned i = 0; i < current_ref_block; ++i) {
			it.push_back(buf[i].begin());
			Join_record::push_next(i, std::numeric_limits<unsigned>::max(), it.back(), records);
		}
		std::make_heap(records.begin(), records.end());
	}
	Target* get_target()
	{
		if (records.empty())
			return 0;
		const unsigned block = records.front().block_;
		const unsigned subject = records.front().info_.subject_id;
		Target *t = new Target(records.front().info_.score);
		do {
			std::pop_heap(records.begin(), records.end());
			records.pop_back();
			if (Join_record::push_next(block, subject, it[block], records))
				std::push_heap(records.begin(), records.end());
		} while (!records.empty() && records.front().info_.subject_id == subject);
	}
	vector<Join_record> records;
	vector<Binary_buffer::Iterator> it;
};

void join_query(vector<Binary_buffer> &buf, Text_buffer &out, Statistics &statistics, unsigned query, const char *query_name, unsigned query_source_len, Output_format &f)
{
	sequence context[6];
	for (unsigned i = 0; i < align_mode.query_contexts; ++i)
		context[i] = query_seqs::get()[query*align_mode.query_contexts + i];

	vector<Join_record> records;
	vector<Binary_buffer::Iterator> it;
	for (unsigned i = 0; i < current_ref_block; ++i) {
		it.push_back(buf[i].begin());
		Join_record::push_next(i, std::numeric_limits<unsigned>::max(), it.back(), records);
	}
	std::make_heap(records.begin(), records.end());
	unsigned block, subject, n_target_seq = 0, hsp_num = 0;
	block = subject = std::numeric_limits<unsigned>::max();
	int top_score = records.front().info_.score;

	while (!records.empty()) {
		const Join_record &next = records.front();
		const unsigned b = next.block_;
		const bool same_subject = n_target_seq > 0 && b == block && next.info_.subject_id == subject;

		if (config.output_range(n_target_seq, next.info_.score, top_score) || same_subject) {
			//printf("q=%u s=%u n=%u ss=%u\n",query, next.info_.subject_id, n_target_seq, same_subject, next.info_.score);

			if(f == Output_format::daa)
				write_daa_record(out, next.info_);
			else {
				Hsp_data hsp(next.info_, query_source_len);
				f.print_match(Hsp_context(hsp,
					query,
					context[align_mode.check_context(hsp.frame)],
					align_mode.query_translated ? query_source_seqs::get()[query] : context[0],
					query_name,
					ref_map.check_id(next.info_.subject_id),
					ref_map.original_id(next.info_.subject_id),
					ref_map.name(next.info_.subject_id),
					ref_map.length(next.info_.subject_id),
					n_target_seq,
					hsp_num).parse().set_query_source_range(next.info_.query_begin), out);
			}

			statistics.inc(Statistics::MATCHES);
			if (!same_subject) {
				block = b;
				subject = next.info_.subject_id;
				++n_target_seq;
				hsp_num = 0;
				statistics.inc(Statistics::PAIRWISE);
			}
			else
				++hsp_num;
		}
		else
			break;

		std::pop_heap(records.begin(), records.end());
		records.pop_back();

		if (Join_record::push_next(block, subject, it[block], records))
			std::push_heap(records.begin(), records.end());
	}
}

void join_worker(Task_queue<Text_buffer,Join_writer> *queue)
{
	Join_fetcher fetcher;
	size_t n;
	Text_buffer *out;
	Statistics stat;
	const String_set<0>& qids = query_ids::get();

	while (queue->get(n, out, fetcher) && fetcher.query_id != Intermediate_record::finished) {
		stat.inc(Statistics::ALIGNED);
		size_t seek_pos;

		const char * query_name = qids[qids.check_idx(fetcher.query_id)].c_str();
		const sequence query_seq = align_mode.query_translated ? query_source_seqs::get()[fetcher.query_id] : query_seqs::get()[fetcher.query_id];

		if (*output_format != Output_format::daa && config.report_unaligned != 0) {
			for (unsigned i = fetcher.unaligned_from; i < fetcher.query_id; ++i) {
				output_format->print_query_intro(i, query_ids::get()[i].c_str(), get_source_query_len(i), *out, true);
				output_format->print_query_epilog(*out, query_ids::get()[i].c_str(), true);
			}
		}

		auto_ptr<Output_format> f(output_format->clone());

		if (*f == Output_format::daa)
			seek_pos = write_daa_query_record(*out, query_name, query_seq);
		else
			f->print_query_intro(fetcher.query_id, query_name, (unsigned)query_seq.length(), *out, false);

		join_query(fetcher.buf, *out, stat, fetcher.query_id, query_name, (unsigned)query_seq.length(), *f);

		if (*f == Output_format::daa)
			finish_daa_query_record(*out, seek_pos);
		else
			f->print_query_epilog(*out, query_name, false);

		queue->push(n);
	}

	statistics += stat;
}

void join_blocks(unsigned ref_blocks, Output_stream &master_out, const vector<Temp_file> &tmp_file)
{
	ref_map.init_rev_map();
	Join_fetcher::init(tmp_file);
	Join_writer writer(master_out);
	Task_queue<Text_buffer, Join_writer> queue(3 * config.threads_, writer);
	Thread_pool threads;
	for (unsigned i = 0; i < config.threads_; ++i)
		threads.push_back(launch_thread(join_worker, &queue));
	threads.join_all();
	Join_fetcher::finish();
	if (*output_format != Output_format::daa && config.report_unaligned != 0) {
		Text_buffer out;
		for (unsigned i = Join_fetcher::query_last + 1; i < query_ids::get().get_length(); ++i) {
			output_format->print_query_intro(i, query_ids::get()[i].c_str(), get_source_query_len(i), out, true);
			output_format->print_query_epilog(out, query_ids::get()[i].c_str(), true);
		}
		writer(out);
	}
}