/****
Copyright (c) 2016, University of Tuebingen, Benjamin Buchfink
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

#ifndef ALIGN_READ_H_
#define ALIGN_READ_H_

#include <vector>
#include <assert.h>
#include "../util/async_buffer.h"
#include "../basic/match.h"
#include "../basic/statistics.h"
#include "align.h"
#include "../util/text_buffer.h"
#include "link_segments.h"
#include "../output/daa_write.h"
#include "../output/output.h"
#include "../output//output_format.h"

// #define ENABLE_LOGGING_AR

#ifdef ENABLE_LOGGING_AR
#include <set>
#endif

using std::vector;

void align_read(Text_buffer &buffer,
		Statistics &stat,
		Trace_pt_buffer::Vector::iterator &begin,
		Trace_pt_buffer::Vector::iterator &end)
{
	static TLS_PTR vector<local_match> *local_ptr = 0;
	static TLS_PTR vector<Segment> *matches_ptr = 0;

#ifdef ENABLE_LOGGING_AR
	static std::set<unsigned> q;
	static tthread::mutex mtx;
#endif

	vector<Segment>& matches (TLS::get(matches_ptr));
	vector<local_match>& local (TLS::get(local_ptr));
	local.clear();
	matches.clear();

	assert(end > begin);
	const size_t hit_count = end - begin;
	local.reserve(hit_count);
	const unsigned contexts = align_mode.query_contexts;
	const unsigned query = begin->query_/contexts;
	const unsigned query_len ((unsigned)query_seqs::data_->length(query*contexts));
	const unsigned source_query_len = align_mode.query_translated ? (unsigned)query_seqs::data_->reverse_translated_len(query*contexts) : query_len;
	const size_t db_letters = ref_header.letters;
	unsigned padding[6];

#ifdef ENABLE_LOGGING_AR
	mtx.lock();
	q.insert(query);
	cout << tthread::thread::get_current_thread_id() << " query=" << query << " n=" << hit_count << endl;
	mtx.unlock();
#endif
	
	typedef Map<vector<hit>::iterator,hit::Query_id<1> > Map_t;
	Map_t hits (begin, end);
	Map_t::Iterator i = hits.begin();
	while(i.valid()) {
		switch (config.local_align_mode) {
		case 1:
			align_sequence_simple(matches, stat, local, padding, db_letters, source_query_len, i.begin(), i.end());
			break;
		default:
			align_sequence_anchored(matches, stat, local, padding, db_letters, source_query_len, i.begin(), i.end());
		}		
		++i;
	}
	
	if(matches.size() == 0)
		return;

	link_segments(matches);

	std::sort(matches.begin(), matches.end());
	unsigned n_hsp = 0, n_target_seq = 0, hit_hsps = 0;
	vector<Segment>::iterator it = matches.begin();
	const int min_raw_score = score_matrix.rawscore(config.min_bit_score == 0
			? score_matrix.bitscore(config.max_evalue, ref_header.letters, query_len) : config.min_bit_score);
	const int top_score = matches.operator[](0).score_;
	size_t seek_pos = 0;

	while(it < matches.end()) {
		const bool same_subject = it != matches.begin() && (it-1)->subject_id_ == it->subject_id_;
		if(!same_subject && it->score_ < min_raw_score)
			break;
		if(!same_subject && !config.output_range(n_target_seq, it->score_, top_score))
			break;
		if(config.local_align_mode == 1 && same_subject && (it-1)->score_ == it->score_) {
			++it;
			continue;
		}
		if(static_cast<double>(it->traceback_->identities)*100/it->traceback_->length < config.min_id
				|| (double)it->traceback_->query_source_range.length()*100/(double)source_query_len < config.query_cover) {
			++it;
			continue;
		}
		if(same_subject && config.single_domain) {
			++it;
			continue;
		}
		
		if (blocked_processing) {
			if (n_hsp == 0)
				seek_pos = Intermediate_record::write_query_intro(buffer, query);
			Intermediate_record::write(buffer, *it, query);
		}
		else {
			if (n_hsp == 0) {
				if (*output_format == Output_format::daa)
					seek_pos = write_daa_query_record(buffer, query_ids::get()[query].c_str(), align_mode.query_translated ? query_source_seqs::get()[query] : query_seqs::get()[query]);
				else
					output_format->print_query_intro(query, query_ids::get()[query].c_str(), source_query_len, buffer);
			}
			if (*output_format == Output_format::daa)
				write_daa_record(buffer, *it, query);
			else
				output_format->print_match(Hsp_context(*it->traceback_,
					query,
					query_seqs::get()[query*contexts + it->frame_],
					query_ids::get()[query].c_str(),
					it->subject_id_,
					ref_ids::get()[it->subject_id_].c_str(),
					(unsigned)ref_seqs::get()[it->subject_id_].length(),
					n_target_seq,
					hit_hsps), buffer);
		}		

		++n_hsp;
		if (!same_subject) {
			++n_target_seq;
			hit_hsps = 0;
		}
		else
			++hit_hsps;
		if(config.alignment_traceback && it->traceback_->gap_openings > 0)
			stat.inc(Statistics::GAPPED);
		stat.inc(Statistics::SCORE_TOTAL, it->traceback_->score);
		++it;
	}

	if (n_hsp > 0) {
		if (!blocked_processing) {
			if (*output_format == Output_format::daa)
				finish_daa_query_record(buffer, seek_pos);
			else
				output_format->print_query_epilog(buffer);
		}
		else
			Intermediate_record::finish_query(buffer, seek_pos);
	}

	stat.inc(Statistics::OUT_MATCHES, matches.size());
	if(!blocked_processing) {
		stat.inc(Statistics::MATCHES, n_hsp);
		stat.inc(Statistics::PAIRWISE, n_target_seq);
		if(n_hsp > 0)
			stat.inc(Statistics::ALIGNED);
	}
	
#ifdef ENABLE_LOGGING_AR
	mtx.lock();
	cout << tthread::thread::get_current_thread_id() << " finish query=" << query << endl;
	q.erase(query);
	for (std::set<unsigned>::const_iterator i = q.begin(); i != q.end(); ++i) {
		cout << *i << ' ';
	}
	cout << endl;
	mtx.unlock();
#endif
}

#endif /* ALIGN_READ_H_ */
