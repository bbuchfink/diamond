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

#include "align_range.h"

Trace_pt_buffer* Trace_pt_buffer::instance;
const Reduction Halfbyte_finger_print::reduction("KR E D Q N C G H LM FY VI W P S T A");
TLS_PTR vector<sequence>* hit_filter::subjects_ptr;

void inner_search2(const sorted_list::const_iterator &q,
	const sorted_list::const_iterator &s,
	Statistics &stats,
	unsigned qbegin,
	unsigned qend,
	unsigned sbegin,
	unsigned send)
{
	for (unsigned i = qbegin; i < qend; ++i) {
		const Letter* query = query_seqs::data_->data(q[i]);
		for (unsigned j = sbegin; j < send; ++j) {
			stats.inc(Statistics::SEED_HITS);
			const Letter* subject = ref_seqs::data_->data(s[j]);
			if (fast_match(query, subject) >= config.min_identities)
				stats.inc(Statistics::TENTATIVE_MATCHES1);
		}
	}
}

template<typename _f>
void tile_search(const sorted_list::const_iterator &q,
	const sorted_list::const_iterator &s,
	Statistics &stats,
	unsigned qbegin,
	unsigned qend,
	unsigned sbegin,
	unsigned send,
	_f f,
	unsigned tile_size)
{
	for (unsigned i = qbegin; i < qend; i += tile_size)
		for (unsigned j = sbegin; j < send;j+=tile_size)
			f(q, s, stats, i, std::min(qend,i+tile_size), j, std::min(send,j+tile_size));
}

// const unsigned levels = 2;
const unsigned tile_size[] = { 1024, 128 };

void register_search(vector<Finger_print>::const_iterator q,
	vector<Finger_print>::const_iterator s,
	Statistics &stats)
{
	stats.inc(Statistics::SEED_HITS, 16);
	Finger_print q1 = *(q++), q2 = *(q++), q3=*(q++), q4=*q, s1=*(s++),s2=*(s++),s3=*(s++),s4=*s;
	if(q1.match(s1) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q1.match(s2) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q1.match(s3) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q1.match(s4) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q2.match(s1) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q2.match(s2) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q2.match(s3) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q2.match(s4) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q3.match(s1) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q3.match(s2) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q3.match(s3) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q3.match(s4) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q4.match(s1) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q4.match(s2) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q4.match(s3) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	if (q4.match(s4) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
}

void query_register_search2(vector<Finger_print>::const_iterator q,
	vector<Finger_print>::const_iterator s,
	vector<Finger_print>::const_iterator s_end,
	Statistics &stats)
{
	Finger_print q1 = *(q++), q2 = *(q++), q3 = *(q++), q4 = *(q++), q5=*(q++),q6=*(q++),q7=*(q++),q8=*q;
	for (; s < s_end; ++s) {
		stats.inc(Statistics::SEED_HITS, 8);
		if (q1.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q2.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q3.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q4.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q5.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q6.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q7.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q8.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	}
}

void query_register_search_for_byte(vector<Finger_print>::const_iterator q,
	vector<Finger_print>::const_iterator s,
	vector<Finger_print>::const_iterator s_end,
	Statistics &stats)
{
	Finger_print q1 = *(q++), q2 = *(q++), q3 = *(q++), q4 = *q;
	const vector<Finger_print>::const_iterator end2 = s_end - (s_end - s) % 4;
	for (; s < end2; ) {
		Finger_print s1 = *(s++), s2 = *(s++), s3=*(s++), s4=*(s++);
		stats.inc(Statistics::SEED_HITS, 4*4);
		if (q1.match(s1) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q2.match(s1) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q3.match(s1) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q4.match(s1) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q1.match(s2) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q2.match(s2) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q3.match(s2) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q4.match(s2) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q1.match(s3) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q2.match(s3) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q3.match(s3) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q4.match(s3) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q1.match(s4) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q2.match(s4) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q3.match(s4) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q4.match(s4) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	}
	for (; s < s_end; ++s) {
		if (q1.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q2.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q3.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q4.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	}
}

void query_register_search(vector<Finger_print>::const_iterator q,
	vector<Finger_print>::const_iterator s,
	vector<Finger_print>::const_iterator s_end,
	Statistics &stats)
{
	Finger_print q1 = *(q++), q2 = *(q++), q3 = *(q++), q4 = *(q++), q5=*(q++),q6=*q;
	const vector<Finger_print>::const_iterator end2 = s_end - (s_end - s) % 4;
	for (; s < end2; ) {
		Finger_print s1 = *(s++), s2 = *(s++), s3 = *(s++), s4 = *(s++);
		stats.inc(Statistics::SEED_HITS, 6 * 4);
		if (q1.match(s1) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q2.match(s1) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q3.match(s1) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q4.match(s1) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q5.match(s1) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q6.match(s1) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q1.match(s2) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q2.match(s2) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q3.match(s2) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q4.match(s2) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q5.match(s2) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q6.match(s2) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q1.match(s3) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q2.match(s3) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q3.match(s3) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q4.match(s3) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q5.match(s3) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q6.match(s3) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q1.match(s4) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q2.match(s4) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q3.match(s4) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q4.match(s4) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q5.match(s4) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q6.match(s4) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	}
	for (; s < s_end; ++s) {
		stats.inc(Statistics::SEED_HITS, 6);
		if (q1.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q2.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q3.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q4.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q5.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (q6.match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	}
}

void query_register_search_simple(vector<Finger_print>::const_iterator q,
	vector<Finger_print>::const_iterator s,
	vector<Finger_print>::const_iterator s_end,
	Statistics &stats)
{
	Finger_print qf[4] = { *(q++), *(q++), *(q++), *(q++) }; // , *(q++), *q
	const vector<Finger_print>::const_iterator end2 = s_end - (s_end - s) % 4;
	char counts[4 * 4];
	unsigned n = 0;
	for (; s < end2; ) {
		Finger_print sf[4] = { *(s++), *(s++), *(s++),*(s++) };		
		stats.inc(Statistics::SEED_HITS, 4 * 4);
		counts[0] += qf[0].match(sf[0]);
		counts[1] += qf[0].match(sf[1]);
		counts[2] += qf[0].match(sf[2]);
		counts[3] += qf[0].match(sf[3]);
		counts[4] += qf[1].match(sf[0]);
		counts[5] += qf[1].match(sf[1]);
		counts[6] += qf[1].match(sf[2]);
		counts[7] += qf[1].match(sf[3]);
		counts[8] += qf[2].match(sf[0]);
		counts[9] += qf[2].match(sf[1]);
		counts[10] += qf[2].match(sf[2]);
		counts[11] += qf[2].match(sf[3]);
		counts[12] += qf[3].match(sf[0]);
		counts[13] += qf[3].match(sf[1]);
		counts[14] += qf[3].match(sf[2]);
		counts[15] += qf[3].match(sf[3]);
		/*counts[16] += qf[4].match(sf[0]);
		counts[17] += qf[4].match(sf[1]);
		counts[18] += qf[4].match(sf[2]);
		counts[19] += qf[4].match(sf[3]);
		counts[20] += qf[5].match(sf[0]);
		counts[21] += qf[5].match(sf[1]);
		counts[22] += qf[5].match(sf[2]);
		counts[23] += qf[5].match(sf[3]);*/
		/*if (qf[i].match(sf[0]) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (qf[i].match(sf[1]) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (qf[i].match(sf[2]) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
		if (qf[i].match(sf[3]) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);*/
		__m128i d = _mm_set_epi8(counts[0], counts[1], counts[2], counts[3], counts[4], counts[5], counts[6], counts[7], counts[8], counts[9], counts[10], counts[11], counts[12], counts[13], counts[14], counts[15]);
		n += _mm_movemask_epi8(d);
	}
	for (; s < s_end; ++s) {
		n += qf[0].match(*s);
		n += qf[1].match(*s);
		n += qf[2].match(*s);
		n += qf[3].match(*s);
		/*n += qf[4].match(*s);
		n += qf[5].match(*s);*/
	}
			//if (qf[i].match(*s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1);
	stats.inc(Statistics::TENTATIVE_MATCHES1, n);
}

void inner_search(vector<Finger_print>::const_iterator q,
	vector<Finger_print>::const_iterator q_end,
	vector<Finger_print>::const_iterator s,
	vector<Finger_print>::const_iterator s_end,
	Statistics &stats)
{
	for (; q < q_end; ++q) {
		for (vector<Finger_print>::const_iterator s2 = s; s2 < s_end; ++s2) {
			stats.inc(Statistics::SEED_HITS);
			if (q->match(*s2) >= config.min_identities)
				stats.inc(Statistics::TENTATIVE_MATCHES1);
		}
	}
}

void tiled_search(vector<Finger_print>::const_iterator q,
	vector<Finger_print>::const_iterator q_end,
	vector<Finger_print>::const_iterator s,
	vector<Finger_print>::const_iterator s_end,
	Statistics &stats,
	unsigned level)
{
	switch (level) {
	case 0:
	case 1:
		for (; q < q_end; q += tile_size[level])
			for (vector<Finger_print>::const_iterator s2 = s; s2 < s_end; s2 += tile_size[level])
				tiled_search(q, q + std::min(q_end - q, (ptrdiff_t)tile_size[level]), s2, s2 + std::min(s_end - s2, (ptrdiff_t)tile_size[level]), stats, level+1);
		break;
	case 2:
		for (; q < q_end; q += 6)
			if (q_end - q < 6)
				inner_search(q, q_end, s, s_end, stats);
			else
				query_register_search(q, s, s_end, stats);

	/*case 2:
		if(q_end-q<4 || s_end-s<4)
			inner_search(q, q_end, s, s_end, stats);
		else
			register_search(q, s, stats);*/
	}	
}

template<typename _f>
void tile_search(vector<Finger_print>::const_iterator q,
	vector<Finger_print>::const_iterator q_end,
	vector<Finger_print>::const_iterator s,
	vector<Finger_print>::const_iterator s_end,
	Statistics &stats,
	_f f,
	unsigned tile_size)
{
	for (; q < q_end; q += tile_size)
		for (vector<Finger_print>::const_iterator s2 = s; s2 < s_end; s2 += tile_size)
			//f(q, std::min(q_end, q + tile_size), s, std::min(s_end, s + tile_size), stats);
			f(q, q+std::min(q_end-q, (ptrdiff_t)tile_size), s2, s2+std::min(s_end-s2, (ptrdiff_t)tile_size), stats);
}

void load_fps(const sorted_list::const_iterator &i, vector<Finger_print> &v, const Sequence_set &seqs)
{
	v.clear();
	v.reserve(i.n);
	for (unsigned j = 0; j < i.n && (i.n<=config.hit_cap || i[j] != 0); ++j)
		v.push_back(Finger_print(seqs.data(i[j]), Masked()));
}

void load_fps2(const sorted_list::const_iterator &i, vector<Finger_print> &v, const Sequence_set &seqs)
{
	v.clear();
	v.reserve(i.n);
	for (unsigned j = 0; j < i.n; ++j)
		v.push_back(Finger_print(seqs.data(i[j])));
}

void search_seed(const sorted_list::const_iterator &q,
	const sorted_list::const_iterator &s,
	Statistics &stats,
	Trace_pt_buffer::Iterator &out,
	const unsigned sid)
{
	//cout << q.n << ' ' << s.n << endl;
	/*if (q.n > config.hit_cap)
		return;*/
	static TLS_PTR vector<Finger_print> *vq_ptr, *vs_ptr;
	vector<Finger_print> &vq(get_tls(vq_ptr)), &vs(get_tls(vs_ptr));
	load_fps2(q, vq, *query_seqs::data_);
	load_fps(s, vs, *ref_seqs::data_);
	//stats.inc(Statistics::TENTATIVE_MATCHES2, vq.size()*vs.size());
	//tile_search(q, s, stats, 0, q.n, 0, s.n, inner_search,128);
	tiled_search(vq.begin(), vq.end(), vs.begin(), vs.end(), stats, 0);
}