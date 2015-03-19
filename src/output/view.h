/****
Copyright (c) 2014, University of Tuebingen
Author: Benjamin Buchfink
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

const unsigned view_buf_size = 32;

struct View_writer
{
	View_writer():
		f_ (program_options::output_file + (program_options::compression==1?".gz":""),
				program_options::compression==1,
				program_options::output_file.length()==0 ? Output_stream::stdout_sink : Output_stream::file_sink)
	{ }
	void operator()(Text_buffer &buf)
	{
		f_.write(buf.get_begin(), buf.size());
		buf.clear();
	}
	Output_stream f_;
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
			if(!daa.read_query_buffer(buf[i]))
				return false;
			else
				++n;
		return true;
	}
	Binary_buffer buf[view_buf_size];
	unsigned n;
	DAA_file &daa;
};

template<typename _val>
struct View_context
{
	View_context(DAA_file &daa, View_writer &writer, const Output_format<_val> &format):
		daa (daa),
		writer (writer),
		queue (3*program_options::threads(), writer),
		format (format)
	{ }
	void operator()(unsigned thread_id)
	{
		try {
			size_t n;
			View_fetcher query_buf (daa);
			Text_buffer *buffer = 0;
			while(!exception_state() && queue.get(n, buffer, query_buf)) {
				for(unsigned j=0;j<query_buf.n;++j) {
					DAA_query_record<_val> r (daa, query_buf.buf[j]);
					for(typename DAA_query_record<_val>::Match_iterator i = r.begin(); i.good(); ++i) {
						if(i->frame > 2 && program_options::forwardonly)
							continue;
						format.print_match(*i, *buffer);
					}
				}
				queue.push(n);
			}
		} catch(std::exception &e) {
			exception_state.set(e);
			queue.wake_all();
		}
	}
	DAA_file &daa;
	View_writer &writer;
	Task_queue<Text_buffer,View_writer> queue;
	const Output_format<_val> &format;
};

template<typename _val>
void view(DAA_file &daa)
{
	score_matrix::instance = auto_ptr<score_matrix> (new score_matrix(daa.score_matrix(),
					daa.gap_open_penalty(),
					daa.gap_extension_penalty(),
					daa.match_reward(),
					daa.mismatch_penalty(),
					_val ()));

	log_stream << "Build version = " << daa.diamond_build() << endl;
	log_stream << "DB sequences = " << daa.db_seqs() << endl;
	log_stream << "DB sequences used = " << daa.db_seqs_used() << endl;
	log_stream << "DB letters = " << daa.db_letters() << endl;

	View_writer writer;
	const Output_format<_val>& format (get_output_format<_val>());
	format.print_header(writer.f_);

	View_context<_val> context (daa, writer, format);
	launch_thread_pool(context, program_options::threads());
}

void view()
{
	DAA_file daa (program_options::daa_file);
	if(daa.mode() == blastn)
		view<Nucleotide>(daa);
	else
		view<Amino_acid>(daa);
}

#endif /* VIEW_H_ */
