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
#include <list>
#include <set>
#include "../basic/sequence.h"
#include "../basic/match.h"
#include "../basic/score_matrix.h"
#include "../search/sse_dist.h"
#include "../align/extend_ungapped.h"
#include "../dp/score_profile.h"
#include "../output/output_format.h"

using std::map;
using std::list;
using std::set;

bool disjoint(list<Hsp_traits>::const_iterator begin, list<Hsp_traits>::const_iterator end, const Hsp_traits &t)
{
	for (; begin != end; ++begin)
		if (!begin->disjoint(t) || !begin->collinear(t))
		//if (!begin->rel_disjoint(t))
			return false;
	return true;
}

bool disjoint(list<Hsp_traits>::const_iterator begin, list<Hsp_traits>::const_iterator end, const Diagonal_segment &d)
{
	for (; begin != end; ++begin)
		if (!begin->disjoint(d) || !begin->collinear(d))
		//if (!begin->rel_disjoint(d))
			return false;
	return true;
}

void Diag_graph::clear_edges()
{
	edges.clear();
	for (vector<Diagonal_node>::iterator i = nodes.begin(); i < nodes.end(); ++i)
		i->reset();
}

void Diag_graph::load(vector<Seed_hit>::const_iterator begin, vector<Seed_hit>::const_iterator end)
{
	int d = std::numeric_limits<int>::min(), max_j_end = d;
	for (vector<Seed_hit>::const_iterator i = begin; i < end; ++i) {
		const int d2 = i->diagonal();
		if (d2 != d) {
			d = d2;
			nodes.push_back(i->ungapped);
			max_j_end = nodes.back().subject_end();
		}
		else if (max_j_end < i->ungapped.j) {
			nodes.push_back(i->ungapped);
			max_j_end = std::max(max_j_end, nodes.back().subject_end());
		}
	}
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
	//nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
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

	int get_approximate_link(int d_idx, int e_idx, double space_penalty, int max_i)
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
			prefix_score = diags.nodes[e_idx].prefix_score + gap_score - int(space_penalty*std::max(space - 1, 0)) + d.score;
			path_max = diags.nodes[e_idx].path_max;
			link_score = e.score + d.score + gap_score;
			link_j = d.j;
		}
		
		if (prefix_score > d.score) {
			diags.add_edge(Diag_graph::Edge(prefix_score, std::max(path_max, prefix_score), link_j, d_idx, e_idx, 0));
			if (log)
				cout << "Link n=" << e_idx << " d=" << e.diag() << " i_end=" << e.query_end() << " max_i=" << max_i << " shift=" << shift << " space=" << space << " prefix_score=" << prefix_score << " link_score=" << link_score << endl;
		}
		return prefix_score;
	}

	template<typename _it>
	void forward_pass(_it begin, _it end, bool init, int max_shift, double space_penalty)
	{
		window.clear();

		for (_it it = begin; it != end; ++it) {

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
				//if (d.j - e.subject_end() > max_dist) {
				if (e.prefix_score - int(space_penalty*(std::max(d.j - e.subject_end(), 0))) <= 0) {
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
				//if(get_approximate_link(node, j->second, space_penalty, max_j) > e.prefix_score)
				get_approximate_link(node, j->second, space_penalty, max_j);
				max_j = std::max(max_j, std::min(d.j, e.subject_end()));
				if (e.subject_end() - (d.subject_end() - std::min(e.diag() - d.diag(), 0)) >= reverse_link_min_overhang) {
					if (log)
						cout << "Computing reverse link node=" << j->second << endl;
					get_approximate_link(j->second, node, space_penalty, max_j);
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
				//if (d.j - e.subject_end() > max_dist && j != i) {
				if (e.prefix_score - int(space_penalty*(std::max(d.j - e.subject_end(),0))) <= 0 && j != i) {
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
				//if (get_approximate_link(node, j->second, space_penalty, max_i) > e.prefix_score)
				get_approximate_link(node, j->second, space_penalty, max_i);
				if (e.i < d.i)
					max_i = std::max(max_i, std::min(e.query_end(), d.i));
				if (e.subject_end() - (d.subject_end() - std::min(e.diag() - d.diag(), 0)) >= reverse_link_min_overhang) {
					if (log)
						cout << "Computing reverse link node=" << j->second << endl;
					get_approximate_link(j->second, node, space_penalty, max_i);
				}
				++j;
			}
			i->second = node;

			if (log)
				cout << endl;
		}
	}

	void backtrace(size_t node, int j_end, Hsp_data *out, Hsp_traits &t, set<unsigned> *node_list) const
	{
		const Diagonal_node &d = diags[node];
		const int dd = d.diag();
		t.d_max = std::max(t.d_max, dd);
		t.d_min = std::min(t.d_min, dd);
		if (node_list)
			node_list->insert(node);
		int j;
		vector<Diag_graph::Edge>::const_iterator f = diags.get_edge(node, j_end);
		if (f != diags.edges.end()) {
			const Diagonal_node &e = diags[f->node_out];
			const int shift = d.diag() - e.diag();
			j = f->j;

			backtrace(f->node_out, shift > 0 ? j : j + shift, out, t, node_list);

			if (out) {
				if (shift > 0) {
					out->transcript.push_back(op_insertion, (unsigned)shift);
					out->length += shift;
				}
				else if (shift < 0) {
					for (int j2 = j + shift; j2 < j; ++j2) {
						out->transcript.push_back(op_deletion, subject[j2]);
						++out->length;
					}
				}
			}
		}
		else {
			if (out) {
				out->query_range.begin_ = d.i;
				out->subject_range.begin_ = d.j;
			}
			t.query_range.begin_ = d.i;
			t.subject_range.begin_ = d.j;
			j = d.j;
		}

		if (out) {
			const int d2 = d.diag();
			if (log) cout << "Backtrace node=" << node << " i=" << d2 + j << "-" << d2 + j_end << " j=" << j << "-" << j_end << endl;
			for (; j < j_end; ++j) {
				const Letter s = subject[j], q = query[d2 + j];
				if (s == q) {
					out->transcript.push_back(op_match);
					++out->identities;
				}
				else
					out->transcript.push_back(op_substitution, s);
				++out->length;
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

	void backtrace(size_t top_node, Hsp_data *out, Hsp_traits &t, set<unsigned> *node_list) const
	{
		Hsp_traits traits;
		if (top_node != Diag_graph::end) {
			if (out) {
				out->transcript.clear();
				out->query_range.end_ = diags[top_node].query_end();
				out->subject_range.end_ = diags[top_node].subject_end();
				out->score = diags.nodes[top_node].prefix_score;
			}
			traits.score = diags.nodes[top_node].prefix_score;
			traits.subject_range.end_ = diags[top_node].subject_end();
			traits.query_range.end_ = diags[top_node].query_end();
			backtrace(top_node, diags[top_node].subject_end(), out, traits, node_list);
		}
		else {
			traits.score = 0;
			if (out)
				out->score = 0;
		}
		if (out)
			out->transcript.push_terminator();
		t = traits;
	}

	int backtrace(list<Hsp_data> &hsps, list<Hsp_traits> &ts, int cutoff, set<unsigned> *node_list) const
	{
		vector<Diagonal_node*> top_nodes;
		for (size_t i = 0; i < diags.nodes.size(); ++i) {
			Diagonal_node &d = diags.nodes[i];
			if (d.prefix_score >= cutoff && d.prefix_score == d.path_max)
				top_nodes.push_back(&d);
		}
		std::sort(top_nodes.begin(), top_nodes.end(), Diagonal_node::cmp_prefix_score);
		int max_score = 0;
		set<unsigned> list2;

		for (vector<Diagonal_node*>::const_iterator i = top_nodes.begin(); i < top_nodes.end(); ++i) {
			Hsp_traits t;
			const size_t node = *i - diags.nodes.data();
			if (disjoint(ts.begin(), ts.end(), **i)) {
				Hsp_data *hsp = 0;
				if (log)
					hsp = new Hsp_data;
				backtrace(node, hsp, t, node_list ? &list2 : 0);
				if (disjoint(ts.begin(), ts.end(), t)) {
					ts.push_back(t);
					if(hsp)
						hsps.push_back(*hsp);
					if (node_list) {
						node_list->insert(list2.begin(), list2.end());
						list2.clear();
					}
					max_score = std::max(max_score, t.score);
				}
				delete hsp;
			}
		}
		return max_score;
	}

	int run(list<Hsp_data> &hsps, list<Hsp_traits> &ts, int band, double space_penalty, int cutoff, bool anschluss)
	{
		diags.sort();
		if (log) {
			diags.print(query, subject);
			cout << endl << endl;
		}

		forward_pass(Index_iterator(0llu), Index_iterator(diags.nodes.size()), true, band, space_penalty);
		auto_ptr<set<unsigned> > node_list(anschluss ? new set<unsigned> : 0);
		int max_score = backtrace(hsps, ts, cutoff, node_list.get());

		if (log && anschluss) {
			for (list<Hsp_data>::iterator i = hsps.begin(); i != hsps.end(); ++i)
				print_hsp(*i, query);
			cout << endl << endl;
		}

		if (anschluss) {
			if (log)
				cout << "Anschluss:" << endl << endl;
			ts.clear();
			hsps.clear();
			diags.clear_edges();
			forward_pass(node_list->begin(), node_list->end(), true, 9999, space_penalty);
			max_score = backtrace(hsps, ts, cutoff, 0);
		}

		if (log) {
			for (list<Hsp_data>::iterator i = hsps.begin(); i != hsps.end(); ++i)
				print_hsp(*i, query);
			cout << endl << "Smith-Waterman:" << endl;
			smith_waterman(query, subject, diags);
			cout << endl << endl;
		}
		return max_score;
	}

	int run(list<Hsp_data> &hsps, list<Hsp_traits> &ts, int band, int cutoff)
	{
		if (ts.empty())
			return 0;
		if(log)
			cout << "***** Scan run n_hsp=" << ts.size() << " cutoff=" << cutoff << endl;
		diags.init();
		ts.sort(Hsp_traits::cmp_diag);
		list<Hsp_traits>::const_iterator i = ts.begin();
		int d_begin = i->d_min - band, d_end = i->d_max + band;
		++i;
		for (; i != ts.end(); ++i) {
			if (i->d_min - band >= d_end) {
				if (log)
					cout << "Scan " << d_begin << '\t' << d_end << '\t' << d_end - d_begin << endl;
				diag_scores.scan_diags(d_begin, d_end, query, subject, qp, query_bc, log, diags.nodes, true);
				d_begin = i->d_min - band;
				d_end = i->d_max + band;
			}
			else
				d_end = std::max(d_end, i->d_max + band);
		}
		if (log)
			cout << "Scan " << d_begin << '\t' << d_end << '\t' << d_end - d_begin << endl;
		diag_scores.scan_diags(d_begin, d_end, query, subject, qp, query_bc, log, diags.nodes, true);
		if (log)
			cout << endl;

		ts.clear();
		hsps.clear();
		return run(hsps, ts, band, config.space_penalty, cutoff, true);
	}

	int run(list<Hsp_data> &hsps, list<Hsp_traits> &ts, vector<Seed_hit>::const_iterator begin, vector<Seed_hit>::const_iterator end, int band)
	{
		if (log)
			cout << "***** Seed hit run " << begin->diagonal() << '\t' << (end - 1)->diagonal() << '\t' << (end - 1)->diagonal() - begin->diagonal() << endl;
		diags.init();
		diags.load(begin, end);
		return run(hsps, ts, band, 0.1, 0, false);
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

int greedy_align(sequence query, const Long_score_profile &qp, const Bias_correction &query_bc, sequence subject, vector<Seed_hit>::const_iterator begin, vector<Seed_hit>::const_iterator end, bool log, list<Hsp_data> &hsps, list<Hsp_traits> &ts)
{
	const int band = config.padding == 0 ? std::min(64, int(query.length()*0.5)) : config.padding;
	Greedy_aligner2 ga(query, qp, query_bc, subject, log);
	return ga.run(hsps, ts, begin, end, band);
}

int greedy_align(sequence query, const Long_score_profile &qp, const Bias_correction &query_bc, sequence subject, bool log, list<Hsp_data> &hsps, list<Hsp_traits> &ts, int cutoff)
{
	const int band = config.padding == 0 ? std::min(64, int(query.length()*0.5)) : config.padding;
	Greedy_aligner2 ga(query, qp, query_bc, subject, log);
	return ga.run(hsps, ts, band, cutoff);
}