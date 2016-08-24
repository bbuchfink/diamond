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

#include <list>
#include <map>
#include "../basic/sequence.h"
#include "../basic/match.h"
#include "../basic/score_matrix.h"
#include "../search/sse_dist.h"
#include "../align/extend_ungapped.h"
#include "../dp/score_profile.h"

using std::list;
using std::map;

struct Intervals
{
	void get_range(int begin, int end, vector<interval> &out)
	{
		for (vector<interval>::iterator i = data_.begin(); i != data_.end(); ++i)
			if(end < (int)i->begin_) {
				out.push_back(interval(begin, end));
				data_.insert(i, interval(begin, end));
				return;
			}
			else if (begin < (int)i->end_) {
				if(begin < (int)i->begin_)
					out.push_back(interval(begin, i->begin_));
				if (end > (int)i->end_)
					out.push_back(interval(i->end_, end));

			}
	}
private:
	vector<interval> data_;
};

struct Link
{
	Link()
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
};

typedef vector<Link> Link_list;

Diagonal_segment score_diagonal(const Letter *query, const Letter *subject)
{
	int i = 0, j = 0, max_score = 0, score = 0, begin = 0, end = 0;
	while (query[i] != '\xff' && subject[i] != '\xff') {
		score += score_matrix(query[i], mask_critical(subject[i]));
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
	return Diagonal_segment(begin, begin, end - begin, max_score);
}

unsigned cmp_16(const Letter *q, const Letter *s)
{
	const __m128i r1 = _mm_loadu_si128((__m128i const*)(q)),
		r2 = _mm_loadu_si128((__m128i const*)(s));
	return _mm_movemask_epi8(_mm_cmpeq_epi8(r1, r2));
}

size_t cmp_count = 0;

void score_diagonal(sequence query, sequence subject, int q, int s, int len, vector<Diagonal_segment> &out)
{
	/*query = (const Letter*)(0xfffffffffffffff0llu & (long long unsigned)query);
	subject = (const Letter*)(0xfffffffffffffff0llu & (long long unsigned)subject);*/
	for (int i = 0; i < len;) {
		++cmp_count;
		if (popcount32(cmp_16(&query[q + i], &subject[s + i])) >= 2) {
			out.push_back(ungapped_extension(s + i + 8, q + i + 8, query, subject));
			i = out.back().query_pos + out.back().len - q;
			//out.push_back(Diagonal_segment());
		}
		i += 16;
	}
}

int score_range(sequence query, sequence subject, int i, int j, int j_end)
{
	int score = 0;
	while (j < j_end) {
		if (i >= (int)query.length()) // || j >= subject.length())
			cout << 'y' << endl;
		score += score_matrix(query[i], subject[j]);
		++i;
		++j;
	}
	return score;
}

void get_ungapped2(sequence query, sequence subject, vector<Diagonal_segment> &out)
{
	/*const int min_score = 15;
	const unsigned ql = query.length(), sl = subject.length();
	for (unsigned i = 0; i < query.length(); ++i) {
		Diagonal_segment d = score_diagonal(&query[i], &subject[0], std::min(ql - i, sl), out);
		if (d.score >= min_score) {
			d.query_pos += i;
			out.push_back(d);
		}
	}
	unsigned i = 1;
	while (subject[i] != '\xff') {
		Diagonal_segment d = score_diagonal(&query[0], &subject[i], std::min(ql, sl - i), out);
		if (d.score >= min_score) {
			d.subject_pos += i;
			out.push_back(d);
		}
		++i;
	}*/
}

void scan_vicinity(sequence q, sequence s, const Diagonal_segment &d, vector<Diagonal_segment> &out)
{
	const int q0 = std::max((int)d.query_pos - 32, 0),
		q1 = std::min((int)d.query_pos + 32, (int)q.length()),
		s0 = std::max((int)d.subject_pos - 32, 0),
		ql = (int)q.length(),
		l = std::min((int)d.subject_pos + (int)d.len + 32, (int)s.length()) - s0;
	for (int i = q0; i < q1; ++i) {
		score_diagonal(q, s, i, s0, std::min(l, ql - i), out);
	}
}

void get_ungapped(sequence query, sequence subject, vector<Diagonal_segment> &d)
{
	typedef map<int, vector<interval> > Dmap;
	Dmap m;
	for (unsigned i = 0; i < d.size(); ++i)
		m[d[i].diag()].push_back(d[i].subject_range());
	for (Dmap::iterator i = m.begin(); i != m.end(); ++i)
		std::sort(i->second.begin(), i->second.end());
	unsigned n = (unsigned)d.size();
	for (unsigned i = 0; i < n; ++i)
		scan_vicinity(query, subject, d[i], d);
	/*n = d.size();
	for (unsigned i = 1; i < n; ++i)
		scan_vicinity(query, subject, d[i], d);*/
}

int get_hgap_link(const Diagonal_segment &d1, const Diagonal_segment &d2, sequence query, sequence subject, Link &l)
{
	const int d = d2.diag() - d1.diag(),
		j2_end = std::min(std::max((int)d2.subject_pos, d1.subject_last() + d + 1) + 1, d2.subject_end());
	int j1 = std::max(std::min(d1.subject_last(), (int)d2.subject_pos - d - 1), (int)d1.subject_pos),
		j2 = j1 + d + 1,
		i1 = d1.query_pos + (j1 - d1.subject_pos),
		i2 = i1 + 1;
	int score1 = score_range(query, subject, d1.query_pos, d1.subject_pos, j1 + 1) - config.gap_open - d*config.gap_extend,
		score2 = score_range(query, subject, i2, j2, d2.subject_pos + d2.len);
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
	if (d1.diag() > d2.diag())
		return get_vgap_link(d1, d2, query, subject, l);
	else
		return get_hgap_link(d1, d2, query, subject, l);
}

void get_links(const vector<Diagonal_segment> &diag, sequence query, sequence subject, vector<Link_list> &links, vector<bool> &is_root, bool log)
{
	for (unsigned i = 0; i < diag.size(); ++i)
		for (unsigned j = 0; j < diag.size();++j) {
			Link l;
			int link_score = get_link(diag[i], diag[j], query, subject, l);
			if (link_score > (int)diag[i].score && link_score > (int)diag[j].score) {
				l.target = j;
				links[i].push_back(l);
				is_root[l.target] = false;
				if (log) cout << i << ' ' << j << ' ' << link_score << endl;
			}
		}
}

void follow_path(unsigned level, unsigned node, vector<Link_list> &links, int score, int subject_pos, sequence query, sequence subject, const vector<Diagonal_segment> &diag, bool log)
{
	if(log) indent(cout, level) << "Node " << node << " Score=" << score << endl;
	const Link_list &l = links[node];
	const Diagonal_segment& d = diag[node];
	const int diff = subject_pos - d.subject_pos;
	if (l.size() == 0) {
		score += score_range(query, subject, d.query_pos + diff, subject_pos, d.subject_pos + d.len);
		if(log) indent(cout, level) << "Final score=" << score << endl << endl;
	}
	for (Link_list::const_iterator i = l.begin(); i != l.end(); ++i) {
		if (i->subject_pos1 < subject_pos)
			continue;
		const int shift = diag[i->target].diag() - d.diag();
		follow_path(level + 1,
			i->target,
			links,
			score + score_range(query, subject, d.query_pos + diff, subject_pos, i->subject_pos1 + 1) - config.gap_open - abs(shift)*config.gap_extend,
			shift < 0 ? i->subject_pos1 + 1 : i->subject_pos1 + shift + 1,
			query,
			subject,
			diag,
			log);
	}
}

void greedy_align2(sequence query, sequence subject, const vector<Diagonal_segment> &sh, bool log)
{
	vector<Diagonal_segment> diag(sh);
	vector<Link_list> links;
	get_ungapped(query, subject, diag);
	if (log) {
		cout << "Diagonals:" << endl;
		for (vector<Diagonal_segment>::const_iterator i = diag.begin(); i != diag.end(); ++i)
			cout << i - diag.begin() << " i=" << i->query_pos << " j=" << i->subject_pos << " score=" << i->score << endl;
		cout << endl;
	}
	//cout << "cmp=" << cmp_count << endl;
	return;

	links.resize(diag.size());
	vector<bool> is_root(diag.size(), true);
	if(log) cout << "Links:" << endl;
	get_links(diag, query, subject, links, is_root, log);
	if(log) cout << endl;

	if(log) cout << "Paths:" << endl;
	for (unsigned i = 0; i < diag.size(); ++i)
		if (is_root[i])
			follow_path(0, i, links, 0, diag[i].subject_pos, query, subject, diag, log);
}

void scan_ahead(sequence q, sequence s, const Diagonal_segment &d, vector<Diagonal_segment> &out)
{
	static const int look_ahead = 64;
	const int q1 = d.query_pos + config.padding,
		l1 = std::min(d.len + look_ahead, (int)s.length() - d.subject_pos),
		ql = (int)q.length(),
		sl = (int)s.length(),
		s1 = d.subject_pos + config.padding,
		l2 = std::min(d.len + look_ahead, (int)q.length() - d.query_pos);
	for (int i = d.query_pos + 1; i < q1; ++i)
		score_diagonal(q, s, i, d.subject_pos, std::min(l1, ql - i), out);
	for (int j = d.subject_pos + 1; j < s1; ++j)
		score_diagonal(q, s, d.query_pos, j, std::min(l2, sl - j), out);
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

struct Diag_buffer
{
	struct Fragment
	{
		Fragment():
			j0(0),
			j1(0),
			score(0)
		{}
		int j0, j1, score;
	};

	void init(int query_pos, int subject_pos)
	{
		this->query_pos = query_pos;
		this->subject_pos = subject_pos;
	}
	uint8_t* get_col(int j)
	{
		if ((int)data_.size() < (j + 1) * 64)
			data_.resize((j + 1) * 64);
		return &data_[j * 64];
	}
	pair<int, int> get_diag(int d, uint8_t max) const
	{
		vector<uint8_t>::const_iterator i = data_.begin() + d, j = i - 64;
		for (; i < data_.end(); i += 64) {
			if (*i == 0)
				j = i;
			if (*i == max)
				break;
		}
		return pair<int, int>((int)(j - data_.begin()) / 64 + 1, (int)(i - data_.begin()) / 64 + 1);
	}
	int query_pos, subject_pos;
	vector<uint8_t> data_;
};

void greedy_align(sequence query, sequence subject, const Diagonal_segment &sh, bool log)
{
	vector<Diagonal_segment> diag;
	Diagonal_segment d = sh;
	while (true) {
		scan_ahead(query, subject, d, diag);
		if (log) {
			cout << "Diagonals:" << endl;
			for (vector<Diagonal_segment>::const_iterator i = diag.begin(); i != diag.end(); ++i)
				cout << i - diag.begin() << " i=" << i->query_pos << " j=" << i->subject_pos << " d=" << (int)i->subject_pos-(int)i->query_pos << " score=" << i->score << endl;
			cout << endl;
		}
		if (log)
			cout << "Links:" << endl;
		Link link;
		unsigned next = std::numeric_limits<unsigned>::max();
		int max_score = 0;
		for (unsigned i = 0; i < diag.size(); ++i) {
			/*if (abs(d.diag() - diag[i].diag())*config.gap_extend + config.gap_open >= diag[i].score)
				continue;*/
			int link_score;
			if ((link_score = get_link(d, diag[i], query, subject, link)) > (int)d.score) {
				if (log)
					cout << i << ' ' << link_score << endl;
				if (link_score > max_score) {
					max_score = link_score;
					next = i;
				}
			}
		}
		if (max_score == 0)
			break;
		d = diag[next];
		diag.clear();
	}
}

void scan_ahead(const Long_score_profile &qp, sequence s, const Diagonal_segment &d, Diag_buffer &scores, vector<Diagonal_segment> &out)
{
	const int l = std::min((int)d.len + 32, (int)s.length() - (int)d.subject_pos);
	score_vector<uint8_t> vbias(score_matrix.bias());
	score_vector<uint8_t> v[4], max[4];
	int i = (int)d.query_pos - 32;
	for (int j = 0; j < l; ++j) {
		const uint8_t *q = qp.get(s[d.subject_pos + j], i);
		uint8_t* dst = scores.get_col(j);
		v[0] = v[0] + score_vector<uint8_t>(q);
		v[0].unbias(vbias);
		max[0].max(v[0]);
		_mm_storeu_si128((__m128i*)dst, v[0].data_);
		dst += 16;
		q += 16;
		v[1] = v[1] + score_vector<uint8_t>(q);
		v[1].unbias(vbias);
		max[1].max(v[1]);
		_mm_storeu_si128((__m128i*)dst, v[1].data_);
		dst += 16;
		q += 17;
		v[2] = v[2] + score_vector<uint8_t>(q);
		v[2].unbias(vbias);
		max[2].max(v[2]);
		_mm_storeu_si128((__m128i*)dst, v[2].data_);
		dst += 16;
		q += 16;
		v[3] = v[3] + score_vector<uint8_t>(q);
		v[3].unbias(vbias);
		max[3].max(v[3]);
		_mm_storeu_si128((__m128i*)dst, v[3].data_);
		++i;
	}
	for (int b = 0; b < 4; ++b)
		for (int i = 0; i < 16; ++i) {
			uint8_t s;
			//if ((s = max[b].data_.m128i_u8[i]) > 12) {
				pair<int, int> e = scores.get_diag(b * 16 + i, s);
				out.push_back(Diagonal_segment((int)d.query_pos - 32 + b * 16 + i + ((b>1) ? 1 : 0) + e.first, d.subject_pos + e.first, e.second - e.first, s));
			//}
		}
}

void greedy_align(sequence query, const Long_score_profile &qp, sequence subject, const Diagonal_segment &sh, bool log)
{
	static TLS_PTR Diag_buffer *scores_ptr;
	static TLS_PTR vector<Diagonal_segment> *diag_ptr;
	vector<Diagonal_segment> &diag (TLS::get(diag_ptr));
	Diagonal_segment d = sh;
	Diag_buffer &scores (TLS::get(scores_ptr));
	int score = 0, score_last = 0;
	Link link;
	while (true) {
		diag.clear();
		scores.data_.clear();
		scan_ahead(qp, subject, d, scores, diag);
		for (vector<Diagonal_segment>::iterator i = diag.begin(); i != diag.end(); ++i) {
			int l;
			i->score += xdrop_ungapped_right(&query[i->query_pos + i->len], &subject[i->subject_pos + i->len], l);
			i->len += l;
		}
		if (log) {
			cout << "Diagonals:" << endl;
			for (vector<Diagonal_segment>::const_iterator i = diag.begin(); i != diag.end(); ++i) {
				cout << i - diag.begin() << " i=" << i->query_pos << " j=" << i->subject_pos << " d=" << (int)i->subject_pos - (int)i->query_pos << " score=" << i->score << endl;
				cout << sequence(query, i->query_pos, i->query_last()) << endl;
				cout << sequence(subject, i->subject_pos, i->subject_last()) << endl;
			}
			cout << endl;
		}
		if (log)
			cout << "Links:" << endl;

		unsigned next = std::numeric_limits<unsigned>::max();
		int max_score = 0;
		for (unsigned i = 0; i < diag.size(); ++i) {
			if (abs(d.diag() - diag[i].diag())*config.gap_extend + config.gap_open >= (int)diag[i].score)
				continue;
			int link_score;
			Link link2;
			if ((link_score = get_link(d, diag[i], query, subject, link2)) > (int)d.score) {
				if (log) {
					cout << i << ' ' << link_score << "(" << link2.score1 << ", " << link2.score2 << ")" << endl;
					cout << sequence(query, d.query_pos, link2.query_pos1) << endl;
					cout << sequence(subject, d.subject_pos, link2.subject_pos1) << endl;
					cout << sequence(query, link2.query_pos2, diag[i].query_last()) << endl;
					cout << sequence(subject, link2.subject_pos2, diag[i].subject_last()) << endl;
				}
				if (link_score > max_score) {
					max_score = link_score;
					next = i;
					link = link2;
				}
			}
		}
		if (log)
			cout << endl;
		if (max_score == 0) {
			score += score_last;
			if (log) {
				cout << "score = " << score << endl;
				cout << sequence(query, d.query_pos, d.query_last()) << endl;
				cout << sequence(subject, d.subject_pos, d.subject_last()) << endl;
			}
			break;
		}
		
		if ((int)diag[next].query_pos > d.query_last() + 2 && (int)diag[next].subject_pos > d.subject_last() + 2) {
			const sequence q1 = sequence(query, d.query_last() + 1, diag[next].query_pos - 1),
				s1 = sequence(subject, d.subject_last() + 1, diag[next].subject_pos - 1);
			if(log)
				cout << q1 << endl << s1 << endl;
			int nw = needleman_wunsch(q1, s1);
			if (log)
				cout << "nw = " << nw << endl;
			score += d.score + nw;
			score_last = diag[next].score;
			d = diag[next];
		}
		else {
			d = Diagonal_segment(link.query_pos2, link.subject_pos2, ((int)diag[next].query_pos - link.query_pos2) + (int)diag[next].len, link.score2);
			score += link.score1;
			score_last = link.score2;
		}

		if (log)
			cout << endl;
	}
}

struct Tile_map
{
	void init(int ql, int sl)
	{
		const int n_tile_rows = (ql + sl - 1 + 15) / 16;
		scores.init(n_tile_rows, (sl + 15) / 16);
		active.init(n_tile_rows, (sl + 15) / 16, true);
	}
	Matrix<__m128i> scores;
	Matrix<bool> active;
};

void scan_diags(sequence q, const Long_score_profile &qp, sequence s, int abs_diag, int j, int l, Tile_map &tile_map, vector<Diagonal_segment> &out)
{

}

void scan_vicinity(sequence q, const Long_score_profile &qp, sequence s, const Diagonal_segment &d, vector<Diagonal_segment> &out)
{
	static const int space = 32, band = 32;
	const int q0 = std::max((int)d.query_pos - band, 0),
		q1 = std::min((int)d.query_pos + band + 1, (int)q.length()),
		s0 = std::max((int)d.subject_pos - space, 0),
		ql = (int)q.length(),
		l = std::min((int)d.subject_pos + (int)d.len + space, (int)s.length()) - s0;
	for (int i = q0; i < q1; ++i) {
		score_diagonal(q, s, i, s0, std::min(l, ql - i), out);
	}
}

void greedy_align(sequence query, const Long_score_profile &qp, sequence subject, const vector<Diagonal_segment> &sh, bool log)
{
	static TLS_PTR Tile_map *tile_map_ptr;
	static TLS_PTR vector<Diagonal_segment> *diag_ptr;
	vector<Diagonal_segment> &diag(TLS::get(diag_ptr));
	diag.clear();
	//Tile_map &tile_map(TLS::get(tile_map_ptr));

}