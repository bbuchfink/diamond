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

// #define _ITERATOR_DEBUG_LEVEL 0

#include "../basic/sequence.h"
#include "../basic/match.h"
#include "../basic/score_matrix.h"
#include "../search/sse_dist.h"
#include "../align/extend_ungapped.h"
#include "../dp/score_profile.h"

struct Link
{
	Link():
		subject_pos1(-1)
	{}
	Link(unsigned target, int query_pos, int subject_pos, int score1, int score2):
		target(target),
		subject_pos1(subject_pos),
		query_pos1(query_pos),
		score1(score1),
		score2(score2)
	{}
	unsigned target;
	int subject_pos1, query_pos1, subject_pos2, query_pos2, score1, score2;
	Link& transpose()
	{
		std::swap(subject_pos1, query_pos1);
		std::swap(subject_pos2, query_pos2);
		return *this;
	}
	void reset()
	{
		subject_pos1 = -1;
	}
};

Diagonal_segment score_diagonal(const Letter *query, const Letter *subject, int qbegin=0,int jbegin=0)
{
	int i = 0, j = 0, max_score = 0, score = 0, begin = 0, end = 0;
	while (query[i] != '\xff' && subject[i] != '\xff') {
		score += score_matrix(query[i], subject[i]);
		if (score <= 0) {
			score = 0;
			j = i + 1;
		}
		if (score > max_score) {
			max_score = score;
			begin = j;
			end = i + 1;
		}
		++i;
	}
	return Diagonal_segment(qbegin+begin, jbegin+begin, end - begin, max_score);
}

int score_range(sequence query, sequence subject, int i, int j, int j_end)
{
	int score = 0;
	while (j < j_end) {
		score += score_matrix(query[i], subject[j]);
		++i;
		++j;
	}
	return score;
}

int get_hgap_link(const Diagonal_segment &d1, const Diagonal_segment &d2, sequence query, sequence subject, Link &l)
{
	const int d = d1.diag() - d2.diag(),
		j2_end = std::min(std::max((int)d2.j, d1.subject_last() + d + 1) + 1, d2.subject_end());
	int j1 = std::max(std::min(d1.subject_last(), (int)d2.j - d - 1), (int)d1.j),
		j2 = j1 + d + 1,
		i1 = d1.i + (j1 - d1.j),
		i2 = i1 + 1;
	if (j2 > d2.subject_last()) {
		l.reset();
		return std::numeric_limits<int>::min();
	}
	int score1 = score_matrix(query[i1], subject[j1]),
		//score2 = score_range(query, subject, i2, j2, d2.j + d2.len);
		score2 = score_range(query, subject, i2, j2, d2.j) + d2.score;
	int max_score = std::numeric_limits<int>::min();
	while (true) {
		if (score1 + score2 > max_score) {
			max_score = score1 + score2;
			l.query_pos1 = i1;
			l.subject_pos1 = j1;
			l.query_pos2 = i2;
			l.subject_pos2 = j2;
			l.score1 = score1;
			l.score2 = score2;
		}
		score2 -= score_matrix(query[i2], subject[j2]);
		++i1; ++i2; ++j1; ++j2;
		if (j2 >= j2_end)
			break;
		score1 += score_matrix(query[i1], subject[j1]);
	}
	return max_score;
}

int get_vgap_link(const Diagonal_segment &d1, const Diagonal_segment &d2, sequence query, sequence subject, Link &l)
{
	int s = get_hgap_link(d1.transpose(), d2.transpose(), subject, query, l);
	l.transpose();
	return s;
}

int get_link(const Diagonal_segment &d1, const Diagonal_segment &d2, sequence query, sequence subject, Link &l)
{
	if (d1.diag() < d2.diag())
		return get_vgap_link(d1, d2, query, subject, l);
	else
		return get_hgap_link(d1, d2, query, subject, l);
}

int needleman_wunsch(sequence q, sequence s)
{
	static int scores1[100], scores2[100], hgap[100];
	int g = -12;
	int *score1 = scores1, *score2 = scores2;
	score1[0] = 0;
	hgap[0] = std::numeric_limits<int>::min() + 1;
	for (int i = 1; i <= (int)q.length(); ++i) {
		score1[i] = g;
		hgap[i] = std::numeric_limits<int>::min() + 1;
		g--;
	}
	for (int j = 1; j <= (int)s.length(); ++j) {
		g = std::numeric_limits<int>::min() + 1;
		score1[0] = j == 1 ? 0 : (-10 - j);
		for (int i = 1; i <= (int)q.length(); ++i) {
			int sc = score_matrix(q[i - 1], s[j - 1]) + score1[i - 1];
			sc = std::max(sc, hgap[i]);
			sc = std::max(sc, g);
			score2[i] = sc;
			sc -= 12;
			g = std::max(g - 1, sc);
			hgap[i] = std::max(hgap[i] - 1, sc);
		}
		std::swap(score1, score2);
	}
	return score1[q.length()];
}

void set_global_max(score_vector<uint8_t> *max, score_vector<uint8_t> *global_max, uint8_t *&local_max)
{
	global_max[0].max(max[0]);
	max[0].store(local_max);
	local_max += 16;
	global_max[1].max(max[1]);
	max[1].store(local_max);
	local_max += 16;
	global_max[2].max(max[2]);
	max[2].store(local_max);
	local_max += 16;
	global_max[3].max(max[3]);
	max[3].store(local_max);
	local_max += 16;
}

void scan_cols(const Long_score_profile &qp, sequence s, int i, int j, int j_end, uint8_t *sv_max, bool log, uint8_t *buf, uint8_t *local_max)
{
	typedef score_vector<uint8_t> Sv;
	const Sv vbias(score_matrix.bias());
	Sv v[4], max[4], global_max[4];
	int n = 0;
	for (; j < j_end; ++j) {
		const uint8_t *q = qp.get(s[j], i);
		v[0] = v[0] + score_vector<uint8_t>(q);
		v[0] -= vbias;
		max[0].max(v[0]);
		_mm_storeu_si128((__m128i*)buf, v[0].data_);
		q += 16;
		buf += 16;
		v[1] = v[1] + score_vector<uint8_t>(q);
		v[1] -= vbias;
		max[1].max(v[1]);
		_mm_storeu_si128((__m128i*)buf, v[1].data_);
		q += 16;
		buf += 16;
		v[2] = v[2] + score_vector<uint8_t>(q);
		v[2] -= vbias;
		max[2].max(v[2]);
		_mm_storeu_si128((__m128i*)buf, v[2].data_);
		q += 16;
		buf += 16;
		v[3] = v[3] + score_vector<uint8_t>(q);
		v[3] -= vbias;
		max[3].max(v[3]);
		_mm_storeu_si128((__m128i*)buf, v[3].data_);
		buf += 16;
		if ((n & 15) == 0) {
			set_global_max(max, global_max, local_max);
		}
		++i;
		++n;
	}
	set_global_max(max, global_max, local_max);
	global_max[0].store(sv_max);
	global_max[1].store(sv_max + 16);
	global_max[2].store(sv_max + 32);
	global_max[3].store(sv_max + 48);
}

struct Greedy_aligner2
{

	enum { band = 64, n_path = 2 };

	struct Diagonal_node : public Diagonal_segment
	{
		Diagonal_node() :
			Diagonal_segment()
		{
			memset(end_score, 0xff, sizeof(end_score)*sizeof(int));
		}
		Diagonal_node(int query_pos, int subject_pos, int len, int score) :
			Diagonal_segment(query_pos, subject_pos, len, score)
		{
			memset(end_score, 0xff, sizeof(end_score)*sizeof(int));
		}
		Diagonal_node(const Diagonal_segment &d):
			Diagonal_segment(d)
		{
			memset(end_score, 0xff, sizeof(end_score)*sizeof(int));
		}
		bool visited() const
		{
			return end_score[0] != -1;
		}
		int end_score[n_path];
		int link_j;
		unsigned next[n_path];
	};

	int correct_left(unsigned node, int j)
	{
		const Diagonal_segment &d = diags[node];
		const int shift = j - d.j;
		return shift < 0 ? score_range(query, subject, d.i + shift, j, d.j) : -score_range(query, subject, d.i, d.j, j);
	}
	
	int correct_right(unsigned node, int j)
	{
		const Diagonal_segment &d = diags[node];
		const int end = d.subject_end(), shift = j - end;
		return shift < 0 ? -score_range(query, subject, d.query_end() + shift, j, d.subject_end()) : score_range(query, subject, d.query_end(), d.subject_end(), j);
	}

	void get_diag(int i, int j, int max_score, int o)
	{
		const uint8_t *p = score_buf.data() + o, *b = p - band, *p_end = p + score_buf.size();
		if (max_score >= 255 - score_matrix.bias()) {
			const int i0 = std::max(i - j, 0), j0 = std::max(j - i, 0);
			diags.push_back(score_diagonal(&query[i0], &subject[j0], i0, j0));
			return;
		}
		for (; p < p_end; p += band) {
			if (*p == 0)
				b = p;
			if (*p == max_score)
				break;
		}
		const int b2 = int(b - (score_buf.data() + o)) / band + 1;
		diags.push_back(Diagonal_segment(i + b2, j + b2, int(p - b) / band, max_score));
	}
	
	void scan_diags(const Diagonal_segment &diag)
	{
		const int d = diag.diag() - band / 2,
			d1 = d + band - 1,
			i = std::max(0, d1) - band + 1,
			j = i - d,
			j1 = std::min((int)query.length() - d, (int)subject.length());
		uint8_t sv[band], sv_max[band];
		memset(sv, 0, band);
		memset(sv_max, 0, band);
		score_buf.resize(band * (j1 - j));
		scan_cols(qp, subject, i, j, j1, sv_max, log, score_buf.data(), local_max.data());
		for (int o = 0; o < band; ++o)
			if (sv_max[o] > 12) {
				//get_diag(i + o, j, sv_max[o], o);
			}
	}

	int follow_path(unsigned level, unsigned node, int score, int subject_pos)
	{
		static const int max_dist = 32, min_nw_space = 11;
		static const float space_penalty = -0.5;
		if (log) indent(cout, level) << "Node " << node << " Score=" << score << endl;
		Diagonal_node& d = diags[node];
		const int cl_score = d.score + correct_left(node, subject_pos);
		if (d.visited()) {
			//const int s = score + d.end_score + score_range(query, subject, d.i + diff, subject_pos, d.link_j + 1);
			const int s = score + d.end_score[0] + cl_score + correct_right(node, d.link_j + 1);
			if (log)
				indent(cout, level) << "Visited node final_score=" << s << endl;
			return s;
		}
		int max_score = 0, current_score;
		for (unsigned k = node + 1; k < diags.size(); ++k)
			if (diags[k].j - d.subject_last() < max_dist) {
				if (abs(diags[k].i - d.query_last()) >= max_dist)
					continue;
				const int shift = diags[k].diag() - d.diag();
				int gap_score = -config.gap_open - abs(shift)*config.gap_extend;
				const int space = shift > 0 ? diags[k].j - d.subject_last() : diags[k].i - d.query_last();
				if (gap_score + diags[k].score + int(space_penalty*std::max(space, 0)) <= 0)
					continue;
				int score2, subject_pos1, subject_pos2;
				if (space >= min_nw_space && abs(shift) > 1) {
					const sequence q1 = sequence(query, d.query_last() + 1, diags[k].i - 1),
						s1 = sequence(subject, d.subject_last() + 1, diags[k].j - 1);
					if (log)
						cout << q1 << endl << s1 << endl;
					gap_score = needleman_wunsch(q1, s1);
					score2 = diags[k].score;
					subject_pos1 = d.subject_last();
					subject_pos2 = diags[k].j;
				}
				else {
					Link l;
					int link_score = get_link(d, diags[k], query, subject, l);
					score2 = l.score2;
					subject_pos1 = l.subject_pos1;
					subject_pos2 = l.subject_pos2;
				}
				if (log)
					indent(cout, level) << "Link n=" << k << " gap_score=" << gap_score << " shift=" << shift << " space=" << space << endl;
				if (gap_score + score2 > 0) {
					//const int c = score_range(query, subject, d.i + diff, subject_pos, subject_pos1 + 1);
					const int c = cl_score + correct_right(node, subject_pos1 + 1);
					int s = follow_path(level + 1,
						k,
						score + c + gap_score,
						subject_pos2);
					if (s > max_score) {
						max_score = s;
						d.next[0] = k;
						current_score = c;
						d.link_j = subject_pos1;
					}
				}
			}
			else
				break;
		if(max_score == 0) {
			//score += score_range(query, subject, d.i + diff, subject_pos, d.j + d.len);
			score += cl_score;
			if (log) indent(cout, level) << "Final score=" << score << endl << endl;
			d.end_score[0] = 0;
			d.link_j = d.subject_end();
			return score;
		}
		else {
			d.end_score[0] = max_score - score - current_score;
			return max_score;
		}
	}

	void follow_path_approximate()
	{
		static const int max_dist = 32;
		static const float space_penalty = -0.5;
		for (unsigned node = 0; node < diags.size(); ++node) {
			Diagonal_node& d = diags[node];
			if (log) cout << "Node " << node << " Score=" << d.score << endl;
			int max_score = 0;
			for (int k = node - 1; k >= 0; --k) {
				const Diagonal_node &e = diags[k];
				if (d.j - e.subject_last() < max_dist) {
					if (abs(d.i - e.query_last()) >= max_dist)
						continue;
					const int shift = d.diag() - e.diag();
					int gap_score = -config.gap_open - abs(shift)*config.gap_extend;
					const int space = shift > 0 ? d.j - e.subject_last() : d.i - e.query_last();
					const int prefix_score = e.end_score[0] + gap_score + (space >= 1 ? e.score + int(space_penalty*std::max(space - 1, 0)) : e.partial_score(abs(space) + 1));
					if (log)
						cout << "Link n=" << k << " shift=" << shift << " space=" << space << " prefix_score=" << prefix_score << endl;
					if (prefix_score > max_score) {
						max_score = prefix_score;
						d.next[0] = k;
						d.end_score[0] = prefix_score;
					}
				}
				else
					break;
			}
			if (log) {
				cout << "Final score=" << max_score + d.score << endl << endl;
			}
		}
	}

	Greedy_aligner2(const sequence &query, const Long_score_profile &qp, const sequence &subject, const vector<Diagonal_segment> &sh, bool log) :
		query(query),
		subject(subject),
		qp(qp),
		log(log),
		score_buf(TLS::get(score_buf_ptr)),
		local_max(TLS::get(local_max_ptr)),
		diags(TLS::get(diags_ptr))
	{
		diags.clear();
		scan_diags(sh[0]);
		std::sort(diags.begin(), diags.end(), Diagonal_segment::cmp_subject_end);
		if (log)
			for (int k = 0; k < (int)diags.size(); ++k) {
				const Diagonal_segment &d = diags[k];
				cout << "Diag n=" << k << " i=" << d.i << " j=" << d.j << " score=" << d.score << " len=" << d.len << endl;
				cout << sequence(query, d.i, d.query_last()) << endl;
				cout << sequence(subject, d.j, d.subject_last()) << endl;
			}
		if(log) cout << endl;
		//follow_path(0, 0, 0, diags[0].j);
		//follow_path_approximate();
	}

	static TLS_PTR vector<uint8_t> *score_buf_ptr;
	static TLS_PTR vector<uint8_t> *local_max_ptr;
	static TLS_PTR vector<Diagonal_node> *diags_ptr;
	const sequence query, subject;
	const Long_score_profile &qp;
	const bool log;
	vector<uint8_t> &score_buf, &local_max;
	vector<Diagonal_node> &diags;

};

TLS_PTR vector<uint8_t> *Greedy_aligner2::score_buf_ptr;
TLS_PTR vector<uint8_t> *Greedy_aligner2::local_max_ptr;
TLS_PTR vector<Greedy_aligner2::Diagonal_node> *Greedy_aligner2::diags_ptr;

void greedy_align2(sequence query, const Long_score_profile &qp, sequence subject, const vector<Diagonal_segment> &sh, bool log)
{
	Greedy_aligner2(query, qp, subject, sh, log);
}