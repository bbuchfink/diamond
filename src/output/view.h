/****
Copyright (c) 2014-2016, University of Tuebingen, Benjamin Buchfink
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

#ifndef VIEW_H_
#define VIEW_H_

#include "../basic/config.h"
#include "../util/binary_file.h"
#include "../util/text_buffer.h"
#include "daa_file.h"
#include "../util/binary_buffer.h"
#include "output_format.h"
#include "../util/task_queue.h"
#include "../basic/score_matrix.h"
#include "../util/thread.h"

const unsigned view_buf_size = 32;

struct View_writer
{
	View_writer():
		f_ (config.compression == 1
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
	View_fetcher(DAA_file &daa):
		daa (daa)
	{ }
	bool operator()()
	{
		n = 0;
		for(unsigned i=0;i<view_buf_size;++i)
			if (!daa.read_query_buffer(buf[i], query_num)) {
				query_num -= n - 1;
				return false;
			} else
				++n;
		query_num -= n - 1;
		return true;
	}
	Binary_buffer buf[view_buf_size];
	unsigned n;
	size_t query_num;
	DAA_file &daa;
};

void view_query(DAA_query_record &r, Text_buffer &out, const Output_format &format)
{
	format.print_query_intro(r.query_num, r.query_name.c_str(), (unsigned)r.query_len(), out, false);
	for (DAA_query_record::Match_iterator i = r.begin(); i.good(); ++i) {
		if (i->frame > 2 && config.forwardonly)
			continue;
		format.print_match(i->context(), out);
	}
	format.print_query_epilog(out, false);
}

struct View_context
{
	View_context(DAA_file &daa, View_writer &writer, const Output_format &format):
		daa (daa),
		writer (writer),
		queue (3*config.threads_, writer),
		format (format)
	{ }
	void operator()(unsigned thread_id)
	{
		try {
			size_t n;
			View_fetcher query_buf (daa);
			Text_buffer *buffer = 0;
			while(queue.get(n, buffer, query_buf)) {
				for (unsigned j = 0; j < query_buf.n; ++j) {
					DAA_query_record r(daa, query_buf.buf[j], query_buf.query_num + j);
					view_query(r, *buffer, format);
				}
				queue.push(n);
			}
		} catch(std::exception &e) {
			std::cout << e.what() << std::endl;
			std::terminate();
		}
	}
	DAA_file &daa;
	View_writer &writer;
	Task_queue<Text_buffer,View_writer> queue;
	const Output_format &format;
};

void view()
{
	DAA_file daa (config.daa_file);
	score_matrix = Score_matrix("", daa.lambda(), daa.kappa(), daa.gap_open_penalty(), daa.gap_extension_penalty());

	message_stream << "Scoring parameters: " << score_matrix << endl;
	verbose_stream << "Build version = " << daa.diamond_build() << endl;
	message_stream << "DB sequences = " << daa.db_seqs() << endl;
	message_stream << "DB sequences used = " << daa.db_seqs_used() << endl;
	message_stream << "DB letters = " << daa.db_letters() << endl;
	
	task_timer timer("Generating output");
	View_writer writer;
	output_format = auto_ptr<Output_format>(get_output_format());

	Binary_buffer buf;
	size_t query_num;
	daa.read_query_buffer(buf, query_num);
	DAA_query_record r(daa, buf, query_num);
	Text_buffer out;
	view_query(r, out, *output_format);
	
	output_format->print_header(*writer.f_, daa.mode(), daa.score_matrix(), daa.gap_open_penalty(), daa.gap_extension_penalty(), daa.evalue(), r.query_name.c_str(), (unsigned)r.query_len());
	writer(out);

	View_context context(daa, writer, *output_format);
	launch_thread_pool(context, config.threads_);
	output_format->print_footer(*writer.f_);
}

#endif /* VIEW_H_ */
