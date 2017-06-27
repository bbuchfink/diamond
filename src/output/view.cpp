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

#include "../basic/config.h"
#include "../util/binary_file.h"
#include "../util/text_buffer.h"
#include "daa_file.h"
#include "../util/binary_buffer.h"
#include "output_format.h"
#include "../util/task_queue.h"
#include "../basic/score_matrix.h"
#include "../util/thread.h"
#include "../data/taxonomy.h"

const unsigned view_buf_size = 32;

struct View_writer
{
	View_writer() :
		f_(config.compression == 1
			? new Compressed_ostream(config.output_file)
			: new Output_stream(config.output_file))
	{ }
	void operator()(Text_buffer &buf)
	{
		f_->write(buf.get_begin(), buf.size());
		buf.clear();
	}
	~View_writer()
	{
		f_->close();
	}
	auto_ptr<Output_stream> f_;
};

struct View_fetcher
{
	View_fetcher(DAA_file &daa) :
		daa(daa)
	{ }
	bool operator()()
	{
		n = 0;
		for (unsigned i = 0; i<view_buf_size; ++i)
			if (!daa.read_query_buffer(buf[i], query_num)) {
				query_num -= n - 1;
				return false;
			}
			else
				++n;
		query_num -= n - 1;
		return true;
	}
	Binary_buffer buf[view_buf_size];
	unsigned n;
	size_t query_num;
	DAA_file &daa;
};

void view_query(DAA_query_record &r, Text_buffer &out, Output_format &format)
{
	auto_ptr<Output_format> f(format.clone());
	f->print_query_intro(r.query_num, r.query_name.c_str(), (unsigned)r.query_len(), out, false);
	DAA_query_record::Match_iterator i = r.begin();
	const unsigned top_score = i.good() ? i->score : 0;
	for (; i.good(); ++i) {
		if (i->frame > 2 && config.forwardonly)
			continue;
		if (!config.output_range(i->hit_num, i->score, top_score))
			break;
		f->print_match(i->context(), out);
	}
	f->print_query_epilog(out, r.query_name.c_str(), false);
}

struct View_context
{
	View_context(DAA_file &daa, View_writer &writer, Output_format &format) :
		daa(daa),
		writer(writer),
		queue(3 * config.threads_, writer),
		format(format)
	{ }
	void operator()(unsigned thread_id)
	{
		try {
			size_t n;
			View_fetcher query_buf(daa);
			Text_buffer *buffer = 0;
			while (queue.get(n, buffer, query_buf)) {
				for (unsigned j = 0; j < query_buf.n; ++j) {
					DAA_query_record r(daa, query_buf.buf[j], query_buf.query_num + j);
					view_query(r, *buffer, format);
				}
				queue.push(n);
			}
		}
		catch (std::exception &e) {
			std::cout << e.what() << std::endl;
			std::terminate();
		}
	}
	DAA_file &daa;
	View_writer &writer;
	Task_queue<Text_buffer, View_writer> queue;
	Output_format &format;
};

void view()
{
	task_timer timer("Loading subject IDs");
	DAA_file daa(config.daa_file);
	score_matrix = Score_matrix("", daa.lambda(), daa.kappa(), daa.gap_open_penalty(), daa.gap_extension_penalty());
	timer.finish();

	message_stream << "Scoring parameters: " << score_matrix << endl;
	verbose_stream << "Build version = " << daa.diamond_build() << endl;
	message_stream << "DB sequences = " << daa.db_seqs() << endl;
	message_stream << "DB sequences used = " << daa.db_seqs_used() << endl;
	message_stream << "DB letters = " << daa.db_letters() << endl;

	init_output();
	taxonomy.init();

	timer.go("Generating output");
	View_writer writer;

	Binary_buffer buf;
	size_t query_num;
	if (daa.read_query_buffer(buf, query_num)) {
		DAA_query_record r(daa, buf, query_num);
		Text_buffer out;
		view_query(r, out, *output_format);

		output_format->print_header(*writer.f_, daa.mode(), daa.score_matrix(), daa.gap_open_penalty(), daa.gap_extension_penalty(), daa.evalue(), r.query_name.c_str(), (unsigned)r.query_len());
		writer(out);

		View_context context(daa, writer, *output_format);
		launch_thread_pool(context, config.threads_);
	}
	else {
		Text_buffer out;
		output_format->print_header(*writer.f_, daa.mode(), daa.score_matrix(), daa.gap_open_penalty(), daa.gap_extension_penalty(), daa.evalue(), "", 0);
		writer(out);
	}

	output_format->print_footer(*writer.f_);
}
