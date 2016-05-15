/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef JOIN_BLOCKS_H_
#define JOIN_BLOCKS_H_

#include <algorithm>
#include <vector>
#include <limits>
#include "output_file.h"

using std::endl;
using std::cout;
using std::vector;

void join_blocks(unsigned ref_blocks, DAA_output &master_out, const vector<Temp_file> &tmp_file)
{
	vector<Block_output*> files;
	vector<Block_output::Iterator> records;
	Block_output::Iterator r;
	for(unsigned i=0;i<ref_blocks;++i) {
		files.push_back(new Block_output (i, tmp_file[i]));
		if(files.back()->next(r, std::numeric_limits<unsigned>::max(), std::numeric_limits<unsigned>::max()))
			records.push_back(r);
	}
	std::make_heap(records.begin(), records.end());
	unsigned query, block, subject, n_target_seq = 0;
	query = block = subject = std::numeric_limits<unsigned>::max();
	int top_score=0;
	Output_buffer buf;
	while(!records.empty()) {
		const Block_output::Iterator &next = records.front();
		const unsigned b = next.block_;

		if(next.info_.query_id != query) {
			if(query != std::numeric_limits<unsigned>::max()) {
				buf.finish_query_record();
				master_out.stream().write(buf.get_begin(), buf.size());
				buf.clear();
			}
			query = next.info_.query_id;
			n_target_seq = 0;
			top_score = next.info_.score;
			statistics.inc(Statistics::ALIGNED);
			buf.write_query_record(query);
		}
		const bool same_subject = n_target_seq > 0 && b == block && next.info_.subject_id == subject;
		if(config.output_range(n_target_seq, next.info_.score, top_score) || same_subject) {
			//printf("q=%u s=%u n=%u ss=%u\n",query, next.info_.subject_id, n_target_seq, same_subject, next.info_.score);
			DAA_output::write_record(buf, next.info_);
			statistics.inc(Statistics::MATCHES);
			if(!same_subject) {
				block = b;
				subject = next.info_.subject_id;
				++n_target_seq;
			}
		} else
			;

		std::pop_heap(records.begin(), records.end());
		records.pop_back();
		if(files[b]->next(r, subject, query)) {
			records.push_back(r);
			std::push_heap(records.begin(), records.end());
		}
	}
	if(query != std::numeric_limits<unsigned>::max()) {
		buf.finish_query_record();
		master_out.stream().write(buf.get_begin(), buf.size());
	}
	for(unsigned i=0;i<ref_blocks;++i) {
		files[i]->close_and_delete();
		delete files[i];
	}
}

#endif /* JOIN_BLOCKS_H_ */
