/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

// #define _ITERATOR_DEBUG_LEVEL 0

#include <map>
#include <list>
#include <set>
#include "../basic/sequence.h"
#include "../basic/match.h"
#include "../basic/score_matrix.h"
//#include "../align/extend_ungapped.h"
#include "../output/output_format.h"
#include "../dp/hsp_traits.h"
#include "chaining.h"

using std::map;
using std::list;
using std::set;

bool disjoint(list<Hsp_traits>::const_iterator begin, list<Hsp_traits>::const_iterator end, const Hsp_traits &t, int cutoff)
{
	for (; begin != end; ++begin)
		if (begin->partial_score(t) < cutoff || !begin->collinear(t))
		//if (!begin->disjoint(t) || !begin->collinear(t))
		//if (!begin->rel_disjoint(t))
			return false;
	return true;
}

bool disjoint(list<Hsp_traits>::const_iterator begin, list<Hsp_traits>::const_iterator end, const Diagonal_segment &d, int cutoff)
{
	for (; begin != end; ++begin)
		if (begin->partial_score(d) < cutoff || !begin->collinear(d))
		//if (!begin->disjoint(d) || !begin->collinear(d))
		//if (!begin->rel_disjoint(d))
			return false;
	return true;
}

void Diag_graph::clear_edges()
{
	edges.clear();
	for (vector<Diagonal_node>::iterator i = nodes.begin(); i < nodes.end(); ++i)
		i->deactivate();
}

void Diag_graph::load(vector<Diagonal_segment>::const_iterator begin, vector<Diagonal_segment>::const_iterator end)
{
	int d = std::numeric_limits<int>::min(), max_j_end = d;
	for (vector<Diagonal_segment>::const_iterator i = begin; i < end; ++i) {
		const int d2 = i->diag();
		if (d2 != d) {
			d = d2;
			nodes.push_back(*i);
			max_j_end = nodes.back().subject_end();
		}
		else if (max_j_end < i->j) {
			nodes.push_back(*i);
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
}

void Diag_graph::prune() {
	vector<Diagonal_node> finished;
	list<Diagonal_node> window;
	for (const Diagonal_node &d : nodes) {
		size_t n = 0;
		for (list<Diagonal_node>::iterator i = window.begin(); i != window.end();) {
			if (i->subject_end() > d.j) {
				if (i->score >= d.score && i->j <= d.j && i->subject_end() >= d.subject_end())
					++n;
				++i;
			}
			else {
				finished.push_back(*i);
				i = window.erase(i);
			}
		}
		if (n <= config.chaining_range_cover)
			window.push_back(d);
	}
	for (const Diagonal_node &d : window)
		finished.push_back(d);
	nodes = std::move(finished);
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
		int prefix_score = 0, link_score = 0, link_j, diff1 = 0, path_max, path_min, prefix_score_begin;
		if (space <= 0) {
			vector<Diag_graph::Edge>::const_iterator edge = diags.get_edge(d_idx, d.j);
			if (edge != diags.edges.end() && edge->prefix_score > e.prefix_score + gap_score + d.score)
				return 0;
			/*if (d.prefix_score > e.prefix_score + gap_score + d.score)
				return 0;*/
			Link link;
			if (get_link(e, d, query, subject, link, link_padding) > 0) {
				diff1 = e.score - link.score1;
				const int prefix_e = diags.prefix_score(e_idx, link.subject_pos1, path_max, path_min);
				prefix_score = prefix_e - diff1 + gap_score + link.score2;
				vector<Diag_graph::Edge>::const_iterator edge = diags.get_edge(d_idx, link.subject_pos2);
				if (edge != diags.edges.end() && edge->prefix_score > prefix_score)
					return 0;
				prefix_score_begin = prefix_score - link.score2;
				path_min = std::min(path_min, prefix_score - link.score2);
				if (prefix_e == path_max) {
					path_max -= diff1;
				}
				link_score = link.score1 + link.score2 + gap_score;
				link_j = link.subject_pos2;
				/*if (log)
					cout << "Link score1=" << link.score1 << " score2=" << link.score2 << " j1=" << link.subject_pos1 << " j2=" << link.subject_pos2 << endl;*/
			}
		}
		else {
			prefix_score = e.prefix_score + gap_score - int(space_penalty*std::max(space - 1, 0)) + d.score;
			vector<Diag_graph::Edge>::const_iterator edge = diags.get_edge(d_idx, d.j);
			if (edge != diags.edges.end() && edge->prefix_score > prefix_score)
				return 0;
			prefix_score_begin = prefix_score - d.score;
			path_max = e.path_max;
			path_min = e.path_min;
			path_min = std::min(path_min, prefix_score - d.score);
			link_score = e.score + d.score + gap_score;
			link_j = d.j;
		}
		
		if (prefix_score > d.score) {
			path_max = std::max(path_max, prefix_score);
			diags.add_edge(Diag_graph::Edge(prefix_score, path_max, link_j, d_idx, e_idx, prefix_score == path_max ? prefix_score : path_min, prefix_score_begin));
			if (log)
				cout << "Link n=" << e_idx << " d=" << e.diag() << " i_end=" << e.query_end() << " max_i=" << max_i << " shift=" << shift << " space=" << space << " prefix_score=" << prefix_score << " link_score=" << link_score << " path_min="<<path_min<<endl;
		}
		return prefix_score;
	}

	template<typename _it>
	void forward_pass(_it begin, _it end, bool init, double space_penalty)
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
				if (e.subject_end() < max_j)
					continue;
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
			while (j != window.end()) {
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
				cout << "Prefix_score=" << d.prefix_score << " path_max=" << d.path_max << " path_min=" << d.path_min << endl << endl;
		}
	}

	bool backtrace(size_t node, int j_end, Hsp *out, Hsp_traits &t, int score_max, int score_min, int max_shift, unsigned &next) const
	{
		const Diagonal_node &d = diags[node];
		vector<Diag_graph::Edge>::const_iterator f = diags.get_edge(node, j_end);
		bool at_end = f >= diags.edges.end();
		const int prefix_score = at_end ? d.score : f->prefix_score;
		if (prefix_score > score_max)
			return false;
		
		int j;
		score_min = std::min(score_min, at_end ? 0 : f->prefix_score_begin);
		
		//if (f != diags.edges.end() && (!stop_at_min || f->path_min == diags[f->node_out].path_min)) {
		if (!at_end) {
			const Diagonal_node &e = diags[f->node_out];
			const int shift = d.diag() - e.diag();
			j = f->j;
			
			if (abs(shift) <= max_shift) {
				const bool bt = backtrace(f->node_out, shift > 0 ? j : j + shift, out, t, score_max, score_min, max_shift, next);
				if (!bt) {
					if (f->prefix_score_begin > score_min)
						return false;
					else
						at_end = true;
				}
			}
			else {
				next = f->node_out;
				at_end = true;
			}
		}

		if (at_end) {
			if (out) {
				out->query_range.begin_ = d.i;
				out->subject_range.begin_ = d.j;
				out->score = score_max - score_min;
			}
			t.query_range.begin_ = d.i;
			t.subject_range.begin_ = d.j;
			t.score = score_max - score_min;
			j = d.j;
		}
		else {
			const Diagonal_node &e = diags[f->node_out];
			const int shift = d.diag() - e.diag();
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

		const int dd = d.diag();
		t.d_max = std::max(t.d_max, dd);
		t.d_min = std::min(t.d_min, dd);
		
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
		return true;
	}

	void backtrace(size_t top_node, Hsp *out, Hsp_traits &t, int max_shift, unsigned &next, int max_j) const
	{
		Hsp_traits traits(frame);
		if (top_node != Diag_graph::end) {
			const Diagonal_node &d = diags[top_node];
			if (out) {
				out->transcript.clear();
				out->query_range.end_ = d.query_end();
				out->subject_range.end_ = d.subject_end();
			}
			traits.subject_range.end_ = d.subject_end();
			traits.query_range.end_ = d.query_end();
			int score_min = d.prefix_score;
			backtrace(top_node, std::min(d.subject_end(), max_j), out, traits, d.prefix_score, score_min, max_shift, next);
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

	int backtrace(size_t top_node, list<Hsp> &hsps, list<Hsp_traits> &ts, list<Hsp_traits>::iterator &t_begin, int cutoff, int max_shift) const
	{
		unsigned next;
		int max_score = 0, max_j = (int)subject.length();
		do {
			Hsp *hsp = log ? new Hsp : 0;
			Hsp_traits t(frame);
			next = std::numeric_limits<unsigned>::max();
			backtrace(top_node, hsp, t, max_shift, next, max_j);
			if (t.score > 0)
				max_j = t.subject_range.begin_;
			if (t.score >= cutoff && disjoint(t_begin, ts.end(), t, cutoff)) {
				if (t_begin == ts.end()) {
					ts.push_back(t);
					t_begin = ts.end();
					t_begin--;
				}
				else
					ts.push_back(t);
				if (hsp)
					hsps.push_back(*hsp);
				max_score = std::max(max_score, t.score);
			}
			delete hsp;
			top_node = next;
		} while (next != std::numeric_limits<unsigned>::max());
		return max_score;
	}

	int backtrace(list<Hsp> &hsps, list<Hsp_traits> &ts, int cutoff, int max_shift) const
	{
		vector<Diagonal_node*> top_nodes;
		for (size_t i = 0; i < diags.nodes.size(); ++i) {
			Diagonal_node &d = diags.nodes[i];
			//cout << "node=" << i << " prefix_score=" << d.prefix_score << " path_max=" << d.path_max << " rel_score=" << d.rel_score() << " cutoff=" << cutoff << endl;
			//if (d.prefix_score >= cutoff && (d.prefix_score == d.path_max || d.prefix_score - d.path_min >= cutoff))
			if(d.rel_score() >= cutoff)
				top_nodes.push_back(&d);
		}
		std::sort(top_nodes.begin(), top_nodes.end(), Diagonal_node::cmp_rel_score);
		int max_score = 0;
		list<Hsp_traits>::iterator t_begin = ts.end();

		for (vector<Diagonal_node*>::const_iterator i = top_nodes.begin(); i < top_nodes.end(); ++i) {
			const size_t node = *i - diags.nodes.data();
			if (log)
				cout << "Backtrace candidate node=" << node << endl;
			if (disjoint(t_begin, ts.end(), **i, cutoff)) {
				if (log)
					cout << "Backtrace node=" << node << " prefix_score=" << (*i)->prefix_score << " rel_score=" << (*i)->rel_score() << endl;
				max_score = std::max(max_score, backtrace(node, hsps, ts, t_begin, cutoff, max_shift));
				if (log)
					cout << endl;
			}
		}
		return max_score;
	}

	int run(list<Hsp> &hsps, list<Hsp_traits> &ts, double space_penalty, int cutoff, int max_shift)
	{
		if (config.chaining_maxnodes > 0) {
			std::sort(diags.nodes.begin(), diags.nodes.end(), Diagonal_segment::cmp_score);
			if (diags.nodes.size() > config.chaining_maxnodes)
				diags.nodes.erase(diags.nodes.begin() + config.chaining_maxnodes, diags.nodes.end());
		}
		if (config.chaining_len_cap > 0.0 && diags.nodes.size() > config.chaining_min_nodes) {
			std::sort(diags.nodes.begin(), diags.nodes.end(), Diagonal_segment::cmp_score);
			const double cap = query.length() * config.chaining_len_cap;
			double total_len = 0.0;
			auto it = diags.nodes.begin();
			while (it < diags.nodes.end() && total_len < cap) {
				total_len += it->len;
				++it;
			}
			diags.nodes.erase(std::max(diags.nodes.begin() + config.chaining_min_nodes, it), diags.nodes.end());
		}
		diags.sort();
		diags.prune();
		if (log) {
			diags.print(query, subject);
			cout << endl << endl;
		}

		forward_pass(Index_iterator(0llu), Index_iterator(diags.nodes.size()), true, space_penalty);
		int max_score = backtrace(hsps, ts, cutoff, max_shift);

		if (log) {
			hsps.sort(Hsp::cmp_query_pos);
			for (list<Hsp>::iterator i = hsps.begin(); i != hsps.end(); ++i)
				print_hsp(*i, TranslatedSequence(query));
			cout << endl << "Smith-Waterman:" << endl;
			smith_waterman(query, subject, diags);
			cout << endl << endl;
		}
		return max_score;
	}

	int run(list<Hsp> &hsps, list<Hsp_traits> &ts, vector<Diagonal_segment>::const_iterator begin, vector<Diagonal_segment>::const_iterator end, int band)
	{
		if (log)
			cout << "***** Seed hit run " << begin->diag() << '\t' << (end - 1)->diag() << '\t' << (end - 1)->diag() - begin->diag() << endl;
		diags.init();
		diags.load(begin, end);
		return run(hsps, ts, 0.1, 19, band);
	}

	Greedy_aligner2(const sequence &query, const sequence &subject, bool log, unsigned frame) :
		query(query),
		subject(subject),
		//query_bc(query_bc),
		log(log),
		frame(frame)
	{
	}

	const sequence query, subject;
	//const Bias_correction &query_bc;
	const bool log;
	const unsigned frame;
	static thread_local Diag_graph diags;
	static thread_local map<int, unsigned> window;

};

thread_local Diag_graph Greedy_aligner2::diags;
thread_local map<int, unsigned> Greedy_aligner2::window;

std::pair<int, list<Hsp_traits>> greedy_align(sequence query, sequence subject, vector<Diagonal_segment>::const_iterator begin, vector<Diagonal_segment>::const_iterator end, bool log, unsigned frame)
{
	const int band = config.chaining_maxgap;
	if (end - begin == 1)
		return { begin->score, { { begin->diag(), begin->diag(), begin->score, (int)frame, begin->query_range(), begin->subject_range() } } };
	Greedy_aligner2 ga(query, subject, log, frame);
	list<Hsp> hsps;
	list<Hsp_traits> ts;
	int score = ga.run(hsps, ts, begin, end, band);
	return std::make_pair(score, std::move(ts));
}