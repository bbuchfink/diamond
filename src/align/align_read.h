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

#ifndef ALIGN_READ_H_
#define ALIGN_READ_H_

#include <vector>
#include <assert.h>
#include "../util/async_buffer.h"
#include "../basic/match.h"
#include "../basic/statistics.h"
#include "../search/align_ungapped.h"
#include "align_sequence.h"
#include "../util/text_buffer.h"
#include "../output/output_buffer.h"

using std::vector;

template<typename _val, typename _locr, typename _locl>
void align_read(Output_buffer<_val> &buffer,
		Statistics &stat,
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &begin,
		typename Trace_pt_buffer<_locr,_locl>::Vector::iterator &end)
{
	static thread_specific_ptr<vector<local_match<_val> > > local_ptr;
	static thread_specific_ptr<vector<Segment<_val> > > matches_ptr;

	Tls<vector<Segment<_val> > > matches (matches_ptr);
	Tls<vector<local_match<_val> > > local (local_ptr);
	local->clear();
	matches->clear();

	assert(end > begin);
	const size_t hit_count = end - begin;
	local->reserve(hit_count);
	const unsigned contexts = query_contexts();
	const unsigned query = begin->query_/contexts;
	const size_t query_len (query_seqs<_val>::data_->length(query*contexts));
	const size_t source_query_len = query_translated() ? query_seqs<_val>::data_->reverse_translated_len(query*contexts) : query_len;
	const size_t db_letters = ref_header.letters;
	unsigned padding[6];

	typedef Map<typename vector<hit<_locr,_locl> >::iterator,typename hit<_locr,_locl>::template Query_id<1> > Map_t;
	Map_t hits (begin, end);
	typename Map_t::Iterator i = hits.begin();
	while(i.valid()) {
		align_sequence<_val,_locr,_locl>(*matches, stat, *local, padding, db_letters, source_query_len, i.begin(), i.end());
		++i;
	}

	if(matches->size() == 0)
		return;

	std::sort(matches->begin(), matches->end());
	unsigned n_hsp = 0, n_target_seq = 0;
	typename vector<Segment<_val> >::iterator it = matches->begin();
	const int min_raw_score = score_matrix::get().rawscore(program_options::min_bit_score == 0
			? score_matrix::get().bitscore(program_options::max_evalue, ref_header.letters, query_len) : program_options::min_bit_score);
	const int top_score = matches->operator[](0).score_;

	while(it < matches->end() && program_options::output_range(n_target_seq, it->score_, top_score) && it->score_ >= min_raw_score) {
		if(it != matches->begin() && (it-1)->subject_id_ == it->subject_id_ && (it-1)->score_ == it->score_) {
			++it;
			continue;
		}
		if(static_cast<double>(it->traceback_->identities_)*100/it->traceback_->len_ < program_options::min_id) {
			++it;
			continue;
		}

		if(n_hsp == 0)
			buffer.write_query_record(query);
		buffer.print_match(*it, source_query_len, query_seqs<_val>::get()[query*contexts + it->frame_], query);
		++n_hsp;

		if(!program_options::long_mode || it == matches->begin() || (it-1)->subject_id_ != it->subject_id_)
			++n_target_seq;
		if(program_options::alignment_traceback && it->traceback_->gap_openings_ > 0)
			stat.inc(Statistics::GAPPED);
		++it;
	}

	if(n_hsp > 0)
		buffer.finish_query_record();

	if(program_options::alignment_traceback)
		for(typename vector<Segment<_val> >::iterator it = matches->begin(); it != matches->end(); ++it)
			if(it->traceback_)
				delete it->traceback_->transcript_;

	stat.inc(Statistics::OUT_MATCHES, matches->size());
	if(ref_header.n_blocks == 1) {
		stat.inc(Statistics::MATCHES, n_hsp);
		if(n_hsp > 0)
			stat.inc(Statistics::ALIGNED);
	}
}

#endif /* ALIGN_READ_H_ */
