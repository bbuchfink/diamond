/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#include "output.h"
#include "../util/temp_file.h"
#include "../data/queries.h"
#include "../output/daa_write.h"
#include "output_format.h"

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
		return rhs.same_subject_ ||  (!rhs.same_subject_ && info_.score < rhs.info_.score);
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

void join_query(vector<Binary_buffer> &buf, Text_buffer &out, Statistics &statistics, unsigned query, const char *query_name, unsigned query_source_len)
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

			if(*output_format == Output_format::daa)
				write_daa_record(out, next.info_);
			else {
				Hsp_data hsp(next.info_, query_source_len);
				output_format->print_match(Hsp_context(hsp,
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
				output_format->print_query_epilog(*out, true);
			}
		}

		if (*output_format == Output_format::daa)
			seek_pos = write_daa_query_record(*out, query_name, query_seq);
		else
			output_format->print_query_intro(fetcher.query_id, query_name, (unsigned)query_seq.length(), *out, false);

		join_query(fetcher.buf, *out, stat, fetcher.query_id, query_name, (unsigned)query_seq.length());

		if (*output_format == Output_format::daa)
			finish_daa_query_record(*out, seek_pos);
		else
			output_format->print_query_epilog(*out, false);

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
			output_format->print_query_epilog(out, true);
		}
		writer(out);
	}
}