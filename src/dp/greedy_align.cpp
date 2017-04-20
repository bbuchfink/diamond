/****
Copyright (c) 2016-2017, Benjamin Buchfink
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

#include <map>
#include "../basic/sequence.h"
#include "../basic/match.h"
#include "../basic/score_matrix.h"
#include "../search/sse_dist.h"
#include "../align/extend_ungapped.h"
#include "../dp/score_profile.h"
#include "../output/output_format.h"

using std::map;

void Diag_graph::load(vector<Seed_hit>::const_iterator begin, vector<Seed_hit>::const_iterator end)
{
	for (vector<Seed_hit>::const_iterator i = begin; i < end; ++i)
		nodes.push_back(i->ungapped);
}

void Diag_graph::print(sequence query, sequence subject) const
{
	for (int k = 0; k < (int)nodes.size(); ++k) {
		const Diagonal_segment &d = nodes[k];
		cout << "Diag n=" << k << " i=" << d.i << " j=" << d.j << " d=" << d.diag() << " score=" << d.score << " len=" << d.len << endl;
		cout << sequence(query, d.i, d.query_last()) << endl;
		cout << sequence(subject, d.j, d.subject_last()) << endl;
	}
}

size_t Diag_graph::top_node() const
{
	int top_score = 0, score;
	size_t top_node = end;
	for (size_t k = 0; k < nodes.size(); ++k)
		if ((score = nodes[k].prefix_score) > top_score) {
			top_node = k;
			top_score = score;
		}
	return top_node;
}

void Diag_graph::sort()
{
	std::sort(nodes.begin(), nodes.end(), Diagonal_segment::cmp_subject);
	nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
}

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

struct Greedy_aligner2
{

	enum { link_padding = 10, reverse_link_min_overhang = 10 };

	int get_approximate_link(int d_idx, int e_idx)
	{
		Diagonal_node &d = diags[d_idx];
		Diagonal_node &e = diags[e_idx];
		const int shift = d.diag() - e.diag();
		int gap_score = shift != 0 ? -score_matrix.gap_open() - abs(shift)*score_matrix.gap_extend() : 0;
		const int space = shift > 0 ? d.j - e.subject_last() : d.i - e.query_last();
		int prefix_score = 0, link_score = 0, link_j, diff1 = 0, path_max;
		if (space <= 0) {
			Link link;
			if (get_link(e, d, query, subject, link, link_padding) > 0) {
				diff1 = e.score - link.score1;
				prefix_score = diags.prefix_score(e_idx, link.subject_pos1, path_max) - diff1 + gap_score + link.score2;
				link_score = link.score1 + link.score2 + gap_score;
				link_j = link.subject_pos2;
				/*if (log)
					cout << "Link score1=" << link.score1 << " score2=" << link.score2 << " j1=" << link.subject_pos1 << " j2=" << link.subject_pos2 << endl;*/
			}
		}
		else {
			prefix_score = diags.nodes[e_idx].prefix_score + gap_score - int(config.space_penalty*std::max(space - 1, 0)) + d.score;
			path_max = diags.nodes[e_idx].path_max;
			link_score = e.score + d.score + gap_score;
			link_j = d.j;
		}
		
		if (prefix_score > d.score) {
			diags.add_edge(Diag_graph::Edge(prefix_score, std::max(path_max, prefix_score), link_j, d_idx, e_idx, 0));
			if (log)
				cout << "Link n=" << e_idx << " d=" << e.diag() << " shift=" << shift << " space=" << space << " prefix_score=" << prefix_score << " link_score=" << link_score << endl;
		}
		return prefix_score;
	}

	template<typename _it>
	void forward_pass(_it begin, _it end, bool init, int max_dist, int max_shift)
	{
		window.clear();

		for (_it it = begin; it < end; ++it) {

			unsigned node = (unsigned)(*it);
			if(init) diags.init(node);
			Diagonal_node& d = diags[node];
			const int dd = d.diag();
			if (log) cout << "Node " << node << " Score=" << d.score << endl;
			map<int, unsigned>::iterator i = window.find(dd), j;
			if (i == window.end())
				i = window.insert(std::make_pair(dd, node)).first;

			j = i;
			int max_j = 0;
			if (i == window.begin())
				goto weiter;
			do {
				--j;
				if (dd - j->first > max_shift)
					break;
				Diagonal_node &e = diags[j->second];
				const int de = j->first, shift = dd - de;
				if (d.j - e.subject_end() > max_dist) {
					map<int, unsigned>::iterator k = j;
					if (j == window.begin()) {
						window.erase(j);
						break;
					}
					else {
						++k;
						window.erase(j);
						j = k;
						continue;
					}
				}
				if (e.subject_end() < max_j) {
					continue;
				}
				if(get_approximate_link(node, j->second) > e.prefix_score)
					max_j = std::max(max_j, std::min(d.j, e.subject_end()));
				if (e.subject_end() - (d.subject_end() - std::min(e.diag() - d.diag(), 0)) >= reverse_link_min_overhang) {
					if (log)
						cout << "Computing reverse link node=" << j->second << endl;
					get_approximate_link(j->second, node);
				}
			} while (j != window.begin());

			weiter:
			j = i;
			if (j->second == node)
				++j;
			int max_i = 0;
			while (j != window.end() && j->first - dd < max_shift) {
				Diagonal_node &e = diags[j->second];
				const int de = j->first, shift = dd - de;
				if (d.j - e.subject_end() > max_dist && j != i) {
					map<int, unsigned>::iterator k = j;
					++k;
					window.erase(j);
					j = k;
					continue;
				}
				if (e.query_end() < max_i) {
					++j;
					continue;
				}
				if (get_approximate_link(node, j->second) > e.prefix_score && e.i < d.i)
					max_i = std::max(max_i, std::min(e.query_end(), d.i));
				if (e.subject_end() - (d.subject_end() - std::min(e.diag() - d.diag(), 0)) >= reverse_link_min_overhang) {
					if (log)
						cout << "Computing reverse link node=" << j->second << endl;
					get_approximate_link(j->second, node);
				}
				++j;
			}
			i->second = node;

			if (log)
				cout << endl;
		}
	}

	void backtrace(size_t node, int j_end, Hsp_data &out, bool transcript, Hsp_traits &t) const
	{
		const Diagonal_node &d = diags[node];
		const int dd = d.diag();
		t.d_max = std::max(t.d_max, dd);
		t.d_min = std::min(t.d_min, dd);
		int j;
		vector<Diag_graph::Edge>::const_iterator f = diags.get_edge(node, j_end);
		if (f != diags.edges.end()) {
			const Diagonal_node &e = diags[f->node_out];
			const int shift = d.diag() - e.diag();
			j = f->j;

			backtrace(f->node_out, shift > 0 ? j : j + shift, out, transcript, t);

			if (transcript) {
				if (shift > 0) {
					out.transcript.push_back(op_insertion, (unsigned)shift);
					out.length += shift;
				}
				else if (shift < 0) {
					for (int j2 = j + shift; j2 < j; ++j2) {
						out.transcript.push_back(op_deletion, subject[j2]);
						++out.length;
					}
				}
			}
		}
		else {
			out.query_range.begin_ = d.i;
			out.subject_range.begin_ = d.j;
			j = d.j;
		}

		if (transcript) {
			const int d2 = d.diag();
			if (log) cout << "Backtrace node=" << node << " i=" << d2 + j << "-" << d2 + j_end << " j=" << j << "-" << j_end << endl;
			for (; j < j_end; ++j) {
				const Letter s = subject[j], q = query[d2 + j];
				if (s == q) {
					out.transcript.push_back(op_match);
					++out.identities;
				}
				else
					out.transcript.push_back(op_substitution, s);
				++out.length;
			}
		}
	}

	unsigned backtrace(unsigned node, int j_end) const
	{
		const Diagonal_node &d = diags[node];
		vector<Diag_graph::Edge>::const_iterator f = diags.get_edge(node, j_end);
		if (f != diags.edges.end()) {
			const Diagonal_node &e = diags[f->node_out];
			const int shift = d.diag() - e.diag();
			return backtrace(f->node_out, shift > 0 ? f->j : f->j + shift);
		}
		else
			return node;
	}

	void backtrace(size_t top_node, Hsp_data &out, Hsp_traits &t) const
	{
		out.transcript.clear();
		Hsp_traits traits;
		if (top_node != Diag_graph::end) {
			out.query_range.end_ = diags[top_node].query_end();
			out.subject_range.end_ = diags[top_node].subject_end();			
			backtrace(top_node, diags[top_node].subject_end(), out, log, traits);
			out.score = diags.nodes[top_node].prefix_score;
		}
		else {
			out.score = 0;
		}
		out.transcript.push_terminator();		
		t = traits;
	}

	size_t anschluss()
	{
		if (log)
			cout << "Anschluss:" << endl;
		vector<Diagonal_node*> top_nodes;
		for (size_t i = 0; i < diags.nodes.size(); ++i)
			if (diags.nodes[i].prefix_score == diags.nodes[i].path_max)
				top_nodes.push_back(&diags.nodes[i]);
		std::sort(top_nodes.begin(), top_nodes.end(), Diagonal_node::cmp_prefix_score);

		Interval_set range;
		vector<unsigned> perimeter;
		const int min_score = score_matrix.gap_open() + score_matrix.gap_extend()*config.padding;
		for (vector<Diagonal_node*>::const_iterator i = top_nodes.begin(); i < top_nodes.end(); ++i) {
			if ((*i)->prefix_score >= min_score && range.max_intersect((*i)->subject_range()) == 0) {
				const unsigned node = (unsigned)(*i - diags.nodes.data()),
					begin = backtrace(node, (*i)->subject_end());
				perimeter.push_back(node);
				perimeter.push_back(begin);
				range.push_back(interval(diags.nodes[begin].j, (*i)->subject_end()));
				if (log)
					cout << "Node = " << node << " prefix_score=" << (*i)->prefix_score << ' ' << (**i) << endl;
			}
		}

		if (log)
			cout << endl << endl;

		if (perimeter.size() > 2) {
			std::sort(perimeter.begin(), perimeter.end());
			forward_pass(perimeter.begin(), perimeter.end(), false, std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
			int max = 0;
			size_t max_node = 0;
			for (vector<unsigned>::const_iterator i = perimeter.begin(); i < perimeter.end();++i)
				if (diags.nodes[*i].prefix_score > max) {
					max = diags.nodes[*i].prefix_score;
					max_node = *i;
				}
			return max_node;
		}
		else
			return top_nodes[0] - diags.nodes.data();
	}

	void run(Hsp_data &out, Hsp_traits &t, int band)
	{
		diags.sort();
		if (log) {
			diags.print(query, subject);
			cout << endl << endl;
		}

		forward_pass(Index_iterator(0llu), Index_iterator(diags.nodes.size()), true, band, band * 2);
		//backtrace(anschluss(), out, t);
		backtrace(diags.top_node(), out, t);

		if (log) {
			print_hsp(out, query);
			cout << endl << "Smith-Waterman:" << endl;
			smith_waterman(query, subject, diags);
			cout << endl << endl;
		}
	}

	void run(int d_begin, int d_end, Hsp_data &out, Hsp_traits &t, int band)
	{
		//cout << d_begin << '\t' << d_end << '\t' << d_end-d_begin << endl;
		diags.init();
		diag_scores.scan_diags(d_begin, d_end, query, subject, qp, query_bc, log, diags.nodes, true);
		run(out, t, band);
	}

	void run(Hsp_data &out, Hsp_traits &t, vector<Seed_hit>::const_iterator begin, vector<Seed_hit>::const_iterator end, int band)
	{
		diags.init();
		diags.load(begin, end);
		run(out, t, band);
	}

	Greedy_aligner2(const sequence &query, const Long_score_profile &qp, const Bias_correction &query_bc, const sequence &subject, bool log) :
		query(query),
		subject(subject),
		qp(qp),
		query_bc(query_bc),
		log(log),
		diag_scores(TLS::get(diag_scores_ptr)),
		diags(TLS::get(diags_ptr)),
		window(TLS::get(window_ptr))
	{
	}

	static TLS_PTR Diag_scores *diag_scores_ptr;
	static TLS_PTR Diag_graph *diags_ptr;
	static TLS_PTR map<int, unsigned> *window_ptr;
	const sequence query, subject;
	const Long_score_profile &qp;
	const Bias_correction &query_bc;
	const bool log;	
	Diag_scores &diag_scores;
	Diag_graph &diags;
	map<int, unsigned> &window;

};

TLS_PTR Diag_scores *Greedy_aligner2::diag_scores_ptr;
TLS_PTR Diag_graph *Greedy_aligner2::diags_ptr;
TLS_PTR map<int, unsigned> *Greedy_aligner2::window_ptr;

void greedy_align(sequence query, const Long_score_profile &qp, const Bias_correction &query_bc, sequence subject, int d_begin, int d_end, bool log, Hsp_data &out, Hsp_traits &traits, vector<Seed_hit>::const_iterator begin, vector<Seed_hit>::const_iterator end, int band)
{
	typedef pair<Hsp_data*, Hsp_traits> Hsp_ref;
	Hsp_data hsp;
	Hsp_traits t;
	Greedy_aligner2 ga(query, qp, query_bc, subject, log);
	if (d_end - d_begin > band * 4) {
		if (log)
			cout << "--- Seed hit run" << endl << endl;
		ga.run(hsp, t, begin, end, band);
		d_begin = t.d_min - band;
		d_end = t.d_max + band;
		if (hsp.score > 0)
			ga.run(d_begin, d_end, hsp, t, band);
	}
	else
		ga.run(d_begin, d_end, hsp, t, band);
	if (hsp.score > 0 && out.score > 0) {
		const Hsp_ref h1 = std::make_pair(&out, traits), h2 = std::make_pair(&hsp, t);
		const pair<Hsp_ref, Hsp_ref> h = hsp.subject_range.begin_ > out.subject_range.begin_ ? std::make_pair(h1, h2) : std::make_pair(h2, h1);
		const int d0 = h.first.first->query_range.end_ - h.first.first->subject_range.end_,
			d1 = h.second.first->query_range.begin_ - h.second.first->subject_range.begin_,
			shift = d1 - d0,
			space = shift > 0 ? (int)h.second.first->subject_range.begin_ - (int)h.first.first->subject_range.end_ : (int)h.second.first->query_range.begin_ - (int)h.first.first->query_range.end_,
			s = -abs(shift)*score_matrix.gap_extend() - score_matrix.gap_open() + out.score + hsp.score - int(config.space_penalty*space);
		if (space >= 0 && s > (int)out.score && s > (int)hsp.score) {
			//cout << "Merge " << endl;
			ga.run(std::min(h.first.second.d_min, h.second.second.d_min), std::max(h.first.second.d_max, h.second.second.d_max) + 1, hsp, t, band);
		}
	}
	if (hsp.score > out.score) {
		out = hsp;
		traits = t;
	}
}

double greedy_align(sequence query, const Long_score_profile &qp, const Bias_correction &query_bc, sequence subject, vector<Seed_hit>::const_iterator begin, vector<Seed_hit>::const_iterator end, bool log, Hsp_data &out, Hsp_traits &traits)
{	
	const int band = config.padding == 0 ? std::min(48, (int)query.length() / 2) : config.padding;
	int d_begin = begin->diagonal() - band, d_end = begin->diagonal() + band;
	vector<Seed_hit>::const_iterator z_begin = begin;
	if (log) {
		cout << "Seed hit i=" << begin->query_pos_ << " j=" << begin->subject_pos_ << " d=" << begin->diagonal() << endl;
	}
	for (vector<Seed_hit>::const_iterator i = begin + 1; i < end; ++i) {
		const int d = i->diagonal();
		if(log)
			cout << "Seed hit i=" << i->query_pos_ << " j=" << i->subject_pos_ << " d=" << i->diagonal() << endl;
		if (d - band >= d_end) {
			greedy_align(query, qp, query_bc, subject, d_begin, d_end, log, out, traits, z_begin, i, band);
			d_begin = d - band;
			d_end = d + band;
			z_begin = i;
		}
		else
			d_end = d + band;
	}

	greedy_align(query, qp, query_bc, subject, d_begin, d_end, log, out, traits, z_begin, end, band);

	if (config.use_smith_waterman) {
		int score;
		needleman_wunsch(query, subject, score, Local(), int());
		assert(score >= out.score);
		out.sw_score = score;
		return pow(score - (int)out.score, 2);
	}
	else
		return 0;
}