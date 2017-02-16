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
#include "../output/output_format.h"

struct Link
{
	Link():
		subject_pos1(-1)
	{}
	Link(unsigned target, int query_pos, int subject_pos, int score1, int score2) :
		subject_pos1(subject_pos),
		query_pos1(query_pos),
		score1(score1),
		score2(score2)
	{}
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
		score1 = 0;
		score2 = 0;
	}
};

Diagonal_segment score_diagonal(const Letter *query, const Letter *subject, int qbegin, int jbegin)
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
	return Diagonal_segment(qbegin + begin, jbegin + begin, end - begin, max_score);
}

Diagonal_segment score_diagonal(const Letter *query, const Letter *subject, int n, int qbegin, int jbegin)
{
	int i = 0, j = 0, max_score = 0, score = 0, begin = 0, end = 0;
	while (i < n && query[i] != '\xff' && subject[i] != '\xff') {
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
	return Diagonal_segment(qbegin + begin, jbegin + begin, end - begin, max_score);
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

int get_hgap_link(const Diagonal_segment &d1, const Diagonal_segment &d2, sequence query, sequence subject, Link &l, int padding)
{
	const int d = d1.diag() - d2.diag(),
		j2_end = std::min(std::max((int)d2.j, d1.subject_last() + d + 1 + padding), d2.subject_last());
	int j1;
	bool space;
	if (d1.subject_last() < d2.j - d - 1) {
		j1 = d1.subject_last();
		space = true;
	}
	else {
		j1 = std::max(d2.j - d - 1 - padding, d1.j);
		space = false;
	}
	int j2 = j1 + d + 1,
		i1 = d1.i + (j1 - d1.j),
		i2 = i1 + 1;
	//cout << "j2=" << j2 << " d2.subject_last=" << d2.subject_last() << endl;
	if (j2 > d2.subject_last()) {
		l.reset();
		return std::numeric_limits<int>::min();
	}
	int score1 = 0,
		//score2 = score_range(query, subject, i2, j2, d2.j + d2.len);
		score2 = score_range(query, subject, i2, j2, d2.j) + d2.score - score_range(query, subject, d2.i, d2.j, j2);
	int max_score = std::numeric_limits<int>::min();
	while (true) {
		//cout << "i1=" << i1 << " j1=" << j1 << " i2=" << i2 << " j2=" << j2 << " score1=" << score1 << " score2=" << score2 << " total=" << score1 + score2 << endl;
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
		if (j2 > j2_end)
			break;
		score1 += score_matrix(query[i1], subject[j1]);
	}
	const int j1_end = j2_end - d;
	if (space)
		l.score1 = d1.score + l.score1;
	else
		l.score1 = d1.score - score_range(query, subject, d1.diag() + j1_end, j1_end, d1.subject_end()) + score_range(query, subject, d1.query_end(), d1.subject_end(), j1_end) - score1 + l.score1;
	return max_score;
}

int get_vgap_link(const Diagonal_segment &d1, const Diagonal_segment &d2, sequence query, sequence subject, Link &l, int padding)
{
	int s = get_hgap_link(d1.transpose(), d2.transpose(), subject, query, l, padding);
	l.transpose();
	return s;
}

int get_link(const Diagonal_segment &d1, const Diagonal_segment &d2, sequence query, sequence subject, Link &l, int padding)
{
	if (d1.diag() < d2.diag())
		return get_vgap_link(d1, d2, query, subject, l, padding);
	else
		return get_hgap_link(d1, d2, query, subject, l, padding);
}

void set_global_max(score_vector<uint8_t> *max, score_vector<uint8_t> *global_max, uint8_t *&local_max)
{
#ifdef __SSE2__
	global_max[0].max(max[0]);
	max[0].store(local_max);
	max[0] = score_vector<uint8_t>();
	local_max += 16;
	global_max[1].max(max[1]);
	max[1].store(local_max);
	max[1] = score_vector<uint8_t>();
	local_max += 16;
	global_max[2].max(max[2]);
	max[2].store(local_max);
	max[2] = score_vector<uint8_t>();
	local_max += 16;
	global_max[3].max(max[3]);
	max[3].store(local_max);
	max[3] = score_vector<uint8_t>();
	local_max += 16;
#endif
}

void scan_cols(const Long_score_profile &qp, sequence s, int i, int j, int j_end, uint8_t *sv_max, bool log, uint8_t *buf, uint8_t *local_max, int block_len)
{
#ifdef __SSE2__
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
		//cout << 's' << v[0] << v[1] << v[2] << v[3] << endl;
		if ((n & 15) == 15) {
			//cout << 'l' << max[0] << max[1] << max[2] << max[3] << endl;
			set_global_max(max, global_max, local_max);
		}
		++i;
		++n;
	}
	if(n % block_len != 0)
		set_global_max(max, global_max, local_max);
	global_max[0].store(sv_max);
	global_max[1].store(sv_max + 16);
	global_max[2].store(sv_max + 32);
	global_max[3].store(sv_max + 48);
	//cout << 'g' << global_max[0] << global_max[1] << global_max[2] << global_max[3] << endl;
#endif
}

struct Greedy_aligner2
{

	enum { band = 64, block_len = 16, link_padding = 5, reverse_link_min_overhang = 10 };

	struct Node_ref
	{
		Node_ref():
			node(0),
			score(0)
		{}
		Node_ref(unsigned node, int score) :
			node(node),
			score(score)
		{}
		unsigned node;
		int score;
		bool operator<(const Node_ref &rhs) const
		{
			return score > rhs.score;
		}
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

	void get_diag(int i, int j, int o, int begin, int end, int max_score)
	{
		const uint8_t *p = score_buf.data() + o + begin*band,
			*p_end = p + block_len*band,
			*b0 = p - band,
			*b = b0;
		for (; p < p_end; p += band) {
			if (*p == 0)
				b = p;
			if (*p == max_score)
				break;
		}
		if (b == b0) {
			b0 = score_buf.data() + o;
			for (; b >= b0; b -= band)
				if (*b == 0)
					break;
		}

		p = score_buf.data() + o + end*band;
		p_end = p + block_len*band;
		for (; p < p_end; p += band) {
			if (*p == max_score)
				break;
		}

		const int b2 = int(b - (score_buf.data() + o)) / band + 1;
		diags.push_back(Diagonal_segment(i + b2, j + b2, int(p - b) / band, max_score));
	}

	void get_diag(int i, int j, int o)
	{
		const uint8_t *p = local_max.data() + o, *p_end = local_max.data() + local_max.size();
		int n = 0, begin = -1, best = -1, best_score = -1, second_best_score = -1, second_best = -1, second_best_begin = -1;
		for (; p < p_end; p += band) {
			if (*p >= config.min_diag_raw_score) {
				if (begin == -1) {
					begin = n;
					best = n;
					best_score = *p;
				} 
				else if (*p > best_score) {
					best = n;
					best_score = *p;
				}
				else if (*p > *(p - band)) {
					if (second_best_begin == -1) {
						second_best_begin = n;
						second_best = n;
						second_best_score = *p;
					}
					else {
						if (*p > second_best_score) {
							second_best_score = *p;
							second_best = n;
						}
					}
				}
			}
			else if (begin >= 0) {
				get_diag(i, j, o, begin * 16, best * 16, best_score);
				if (second_best_begin >= 0) {
					second_best_begin -= 1;
					const Diagonal_segment second = score_diagonal(&query[i + second_best_begin * 16],
						&subject[j + second_best_begin * 16],
						(second_best - second_best_begin + 1) * 16,
						i + second_best_begin * 16,
						j + second_best_begin * 16);
					if (second.score >= config.min_diag_raw_score)
						diags.push_back(second);
					second_best_begin = -1;
				}
				begin = -1;
				best = -1;
				best_score = -1;
			}
			++n;
		}
		if (begin >= 0) {
			get_diag(i, j, o, begin * 16, best * 16, best_score);
			if (second_best_begin >= 0) {
				second_best_begin -= 1;
				const Diagonal_segment second = score_diagonal(&query[i + second_best_begin * 16],
					&subject[j + second_best_begin * 16],
					(second_best - second_best_begin + 1) * 16,
					i + second_best_begin * 16,
					j + second_best_begin * 16);
				if (second.score >= config.min_diag_raw_score)
					diags.push_back(second);
			}
		}
	}
	
	void scan_diags(const Diagonal_segment &diag)
	{
		const int d = diag.diag() - band / 2,
			d1 = d + band - 1,
			i = std::max(0, d1) - band + 1,
			j = i - d,
			j1 = std::min((int)query.length() - d, (int)subject.length());
		uint8_t sv_max[band];
		memset(sv_max, 0, band);
		const size_t cells = band * (j1 - j);
		score_buf.clear();
		score_buf.resize(cells);
		memset(score_buf.data(), 0, score_buf.size());
		local_max.clear();
		local_max.resize((j1 - j + block_len - 1) / block_len * band);
		memset(local_max.data(), 0, local_max.size());
		scan_cols(qp, subject, i, j, j1, sv_max, log, score_buf.data(), local_max.data(), block_len);
		for (int o = 0; o < band; ++o)
			if (sv_max[o] >= config.min_diag_raw_score) {
				//get_diag(i + o, j, sv_max[o], o);
				if (sv_max[o] >= 255 - score_matrix.bias()) {
					const int s = std::min(i + o, 0);
					diags.push_back(score_diagonal(&query[i + o - s], &subject[j - s], i + o - s, j - s));
				}
				else
					get_diag(i + o, j, o);
			}
	}

	int follow_path(unsigned level, unsigned node, int score, int subject_pos)
	{
		static const int max_dist = 32, min_nw_space = 7;
		static const float space_penalty = -0.5;
		Diagonal_node& d = diags[node];
		if (log) indent(cout, level) << "Node " << node << " Score=" << d.edges[0].prefix_score << " j=" << subject_pos << endl;
		if (d.diff != std::numeric_limits<int>::min()) {
			if (log)
				indent(cout, level) << "Visited node final_score=" << d.edges[0].prefix_score << endl;
			return d.diff;
		}
		int max_score = d.edges[0].prefix_score, diff = std::numeric_limits<int>::min(), max_edge = 0;
		for (unsigned k = 0; k < Diagonal_node::n_path; ++k) {
			Diagonal_node::Edge &f = d.edges[k];
			if (f.prefix_score == 0)
				break;
			if (f.node == node) {
				if (diff < f.prefix_score - max_score) {
					diff = f.prefix_score - max_score;
					max_edge = k;
				}
				continue;
			}
			if (f.j > subject_pos)
				continue;
			unsigned next = f.node;
			const Diagonal_node &e = diags[next];
			const int shift = d.diag() - e.diag();
			if (!f.exact) {
				int gap_score = -config.gap_open - abs(shift)*config.gap_extend;
				const int space = shift > 0 ? d.j - e.subject_last() : d.i - e.query_last();
				f.prefix_score += int(config.raw_space_penalty*std::max(space - 1, 0));
				if (log)
					indent(cout, level) << "Node=" << node << " Link n=" << f.node << " space=" << space << " shift=" << shift << endl;
				if (space >= min_nw_space && abs(shift) > 1) {
					if (log) {
						const sequence q1 = sequence(query, e.query_last() + 1, d.i - 1),
							s1 = sequence(subject, e.subject_last() + 1, d.j - 1);
						indent(cout, level) << q1 << endl;
						indent(cout, level) << s1 << endl;
					}
					gap_score = needleman_wunsch(query, subject, e.query_last() + 1, d.i, e.subject_last() + 1, d.j, node, k, diags, log);
					f.prefix_score -= -config.gap_open - abs(shift)*config.gap_extend;
					f.prefix_score += gap_score;
				}
				else {
					Link l;
					int link_score = get_link(e, d, query, subject, l, link_padding);
					f.prefix_score += (l.score1 - e.score) + (l.score2 - d.score);
					f.j = l.subject_pos2;
					if (log)
						indent(cout, level) << "Direct link scores=" << l.score1 << "," << l.score2 << endl;
				}
			}
			if (log)
				indent(cout, level) << "Node=" << node << " Link n=" << next << " prefix_score=" << f.prefix_score << " j=" << f.j << endl;
			f.prefix_score += follow_path(level + 1, next, 0, f.j + std::min(shift, 0) - 1);
			if (log)
				indent(cout, level) << "Node=" << node << " Link n=" << next << " final prefix_score=" << f.prefix_score << endl;
			if (diff < f.prefix_score - max_score) {
				diff = f.prefix_score - max_score;
				max_edge = k;
			}
		}
		if (max_edge != 0) {
			memcpy(&d.edges[0], &d.edges[max_edge], sizeof(Diagonal_node::Edge));
		}
		d.diff = diff;
		top_node = std::min(top_node, Node_ref(node, d.edges[0].prefix_score));
		if(log) indent(cout, level) << "Node " << node << " diff=" << diff << " final_score=" << d.edges[0].prefix_score << endl;
		return diff;
	}

	void get_approximate_link(Diagonal_node &d, const Diagonal_node &e, int k)
	{
		const int shift = d.diag() - e.diag();
		int gap_score = -config.gap_open - abs(shift)*config.gap_extend;
		const int space = shift > 0 ? d.j - e.subject_last() : d.i - e.query_last();
		int prefix_score = 0, link_score = 0, link_j;
		bool exact;
		if (space <= 0) {
			Link link;
			if (get_link(e, d, query, subject, link, link_padding) > 0) {
				prefix_score = e.edges[0].prefix_score - (e.score - link.score1) + gap_score + link.score2;
				link_score = link.score1 + link.score2 + gap_score;
				exact = true;
				link_j = link.subject_pos2;
			}
		}
		else {
			prefix_score = e.edges[0].prefix_score + gap_score - int(config.raw_space_penalty*std::max(space - 1, 0)) + d.score;
			link_score = e.score + d.score + gap_score;
			exact = false;
			link_j = d.j;
		}

		if (log)
			cout << "Link n=" << k << " shift=" << shift << " space=" << space << " prefix_score=" << prefix_score << " link_score=" << link_score << endl;
		if (prefix_score > 0)
			d.edges.add(Diagonal_node::Edge(prefix_score, link_j, k, exact));
	}

	unsigned follow_path_approximate()
	{
		static const int max_dist = 60;
		static const float space_penalty = -0.5;
		int max_score = 0;
		unsigned max_node;
		for (unsigned node = 0; node < diags.size(); ++node) {
			Diagonal_node& d = diags[node];
			if (log) cout << "Node " << node << " Score=" << d.score << endl;
			d.edges.add(Diagonal_node::Edge(d.score, d.j, node, true));
			for (int k = node - 1; k >= 0; --k) {
				const Diagonal_node &e = diags[k];
				if (abs(d.j - e.subject_last()) < max_dist) {
					if (abs(d.i - e.query_last()) >= max_dist)
						continue;
					get_approximate_link(d, e, k);
				}
				else
					break;
			}
			if (log) {
				cout << "Final score=" << d.edges[0].prefix_score << endl << endl;
			}
			top_nodes.push_back(Node_ref(node, d.edges[0].prefix_score));
			if (d.edges[0].prefix_score > max_score) {
				max_score = d.edges[0].prefix_score;
				max_node = node;
			}

			for (int k = 0; k < Diagonal_node::n_path; ++k) {
				Diagonal_node::Edge &f = d.edges[k];
				if (f.prefix_score == 0)
					break;
				Diagonal_node &e = diags[f.node];
				if (e.subject_end() - (d.subject_end() - std::min(e.diag()-d.diag(), 0)) >= reverse_link_min_overhang) {
					if (log)
						cout << "Computing reverse link" << endl;
					get_approximate_link(e, d, node);
				}
			}
		}
		return max_node;
	}

	void backtrace(unsigned node, int j_end, Hsp_data &out)
	{
		const Diagonal_node &d = diags[node];
		const Diagonal_node::Edge &f = d.edges[0];
		const Diagonal_node &e = diags[f.node];
		const int shift = d.diag() - e.diag();
		int j = f.j;
		if (f.node != node) {
			backtrace(f.node, shift > 0 ? j : j + shift, out);
			if (shift > 0) {
				out.transcript.push_back(op_insertion, (unsigned)shift);
			}
			else {
				for (int j2 = j + shift; j2 < j; ++j2)
					out.transcript.push_back(op_deletion, subject[j2]);
			}
		}
		else {
			out.query_range.begin_ = d.i;
			out.subject_range.begin_ = d.j;
		}
		const int d2 = d.diag();
		if (log) cout << "Backtrace node=" << node << " i=" << d2+j << "-" << d2+j_end << " j=" << j << "-" << j_end << endl;		
		for (; j < j_end; ++j) {
			const Letter s = subject[j], q = query[d2 + j];
			if (s == q)
				out.transcript.push_back(op_match);
			else
				out.transcript.push_back(op_substitution, s);
		}
	}

	Greedy_aligner2(const sequence &query, const Long_score_profile &qp, const sequence &subject, const vector<Diagonal_segment> &sh, bool log, Hsp_data &out) :
		query(query),
		subject(subject),
		qp(qp),
		log(log),
		score_buf(TLS::get(score_buf_ptr)),
		local_max(TLS::get(local_max_ptr)),
		diags(TLS::get(diags_ptr)),
		top_nodes(TLS::get(top_nodes_ptr))
	{
		static const double score_range = 0.97;

		config.min_diag_raw_score = 15;
		diags.clear();
		local_max.clear();
		top_nodes.clear();
		scan_diags(sh[0]);
		std::sort(diags.begin(), diags.end(), Diagonal_segment::cmp_heuristic);
		if (log)
			for (int k = 0; k < (int)diags.size(); ++k) {
				const Diagonal_segment &d = diags[k];
				cout << "Diag n=" << k << " i=" << d.i << " j=" << d.j << " score=" << d.score << " len=" << d.len << endl;
				cout << sequence(query, d.i, d.query_last()) << endl;
				cout << sequence(subject, d.j, d.subject_last()) << endl;
			}
		if(log) cout << endl << endl;
		unsigned max_node = follow_path_approximate();

		std::sort(top_nodes.begin(), top_nodes.end());
		for (vector<Node_ref>::const_iterator i = top_nodes.begin(); i < top_nodes.end() && (double)i->score / top_nodes.begin()->score >= score_range; ++i)
			follow_path(0, i->node, 0, diags[i->node].subject_last());
		if (log) cout << endl << endl;

		out.transcript.clear();
		out.length = 1;
		backtrace(top_node.node, diags[top_node.node].subject_end(), out);
		out.transcript.push_terminator();
		out.score = top_node.score;

		if (log) {
			Text_buffer buf;
			Pairwise_format().print_match(Hsp_context(out, 0, query, query, "", 0, 0, "", 0, 0, 0), buf);
			buf << '\0';
			cout << buf.get_begin();
		}
	}

	static TLS_PTR vector<uint8_t> *score_buf_ptr;
	static TLS_PTR vector<uint8_t> *local_max_ptr;
	static TLS_PTR vector<Diagonal_node> *diags_ptr;
	static TLS_PTR vector<Node_ref> *top_nodes_ptr;
	const sequence query, subject;
	const Long_score_profile &qp;
	const bool log;
	vector<uint8_t> &score_buf, &local_max;
	vector<Diagonal_node> &diags;
	vector<Node_ref> &top_nodes;
	Node_ref top_node;

};

TLS_PTR vector<uint8_t> *Greedy_aligner2::score_buf_ptr;
TLS_PTR vector<uint8_t> *Greedy_aligner2::local_max_ptr;
TLS_PTR vector<Diagonal_node> *Greedy_aligner2::diags_ptr;
TLS_PTR vector<Greedy_aligner2::Node_ref> *Greedy_aligner2::top_nodes_ptr;

void greedy_align2(sequence query, const Long_score_profile &qp, sequence subject, const vector<Diagonal_segment> &sh, bool log, Hsp_data &out)
{
	Greedy_aligner2(query, qp, subject, sh, log, out);
}