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
#include <iostream>
#include "../basic/sequence.h"
#include "../basic/match.h"
#include "../stats/score_matrix.h"
//#include "../align/extend_ungapped.h"
#include "../output/output_format.h"
#include "../util/hsp/approx_hsp.h"
#include "chaining.h"
#include "../util/util.h"
#include "../dp/ungapped.h"
#include "aligner.h"
#include "diag_graph.h"

using std::endl;
using std::cout;
using std::map;
using std::list;
using std::set;
using std::max;
using std::min;
using std::vector;

static const double SPACE_PENALTY = 0.1;

void DiagGraph::clear_edges()
{
	edges.clear();
	for (vector<DiagonalNode>::iterator i = nodes.begin(); i < nodes.end(); ++i)
		i->deactivate();
}

void DiagGraph::load(vector<DiagonalSegment>::const_iterator begin, vector<DiagonalSegment>::const_iterator end)
{
	int d = std::numeric_limits<int>::min(), max_j_end = d;
	for (vector<DiagonalSegment>::const_iterator i = begin; i < end; ++i) {
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

void DiagGraph::print(Sequence query, Sequence subject) const
{
	for (int k = 0; k < (int)nodes.size(); ++k) {
		const DiagonalSegment &d = nodes[k];
		cout << "Diag n=" << k << " i=" << d.i << " j=" << d.j << " d=" << d.diag() << " score=" << d.score << " len=" << d.len << endl;
		cout << Sequence(query, d.i, d.query_last()) << endl;
		cout << Sequence(subject, d.j, d.subject_last()) << endl;
	}
}

size_t DiagGraph::top_node() const
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

void DiagGraph::sort()
{
	std::sort(nodes.begin(), nodes.end(), DiagonalSegment::cmp_subject);
}

void DiagGraph::prune() {
	vector<DiagonalNode> finished;
	list<DiagonalNode> window;
	for (const DiagonalNode&d : nodes) {
		size_t n = 0;
		for (list<DiagonalNode>::iterator i = window.begin(); i != window.end();) {
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
	for (const DiagonalNode&d : window)
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

int get_hgap_link(const DiagonalSegment &d1, const DiagonalSegment &d2, Sequence query, Sequence subject, Link &l, int padding)
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
		l.score1 += d1.score;
	else
		l.score1 += d1.score - score_range(query, subject, d1.diag() + j1_end, j1_end, d1.subject_end()) + score_range(query, subject, d1.query_end(), d1.subject_end(), j1_end) - score1;
	return max_score;
}

int get_vgap_link(const DiagonalSegment &d1, const DiagonalSegment &d2, Sequence query, Sequence subject, Link &l, int padding)
{
	int s = get_hgap_link(d1.transpose(), d2.transpose(), subject, query, l, padding);
	l.transpose();
	return s;
}

int get_link(const DiagonalSegment &d1, const DiagonalSegment &d2, Sequence query, Sequence subject, Link &l, int padding)
{
	if (d1.diag() < d2.diag())
		return get_vgap_link(d1, d2, query, subject, l, padding);
	else
		return get_hgap_link(d1, d2, query, subject, l, padding);
}

namespace Chaining {

	enum { link_padding = 10, reverse_link_min_overhang = 10 };

	int Aligner::get_approximate_link(int d_idx, int e_idx, double space_penalty, int max_i)
	{
		DiagonalNode& d = diags[d_idx];
		DiagonalNode& e = diags[e_idx];
		const int shift = d.diag() - e.diag();
		int gap_score = shift != 0 ? -score_matrix.gap_open() - abs(shift)*score_matrix.gap_extend() : 0;
		const int space = shift > 0 ? d.j - e.subject_last() : d.i - e.query_last();
		int prefix_score = 0, link_score = 0, link_j, diff1 = 0, path_max, path_min, prefix_score_begin;
		if (space <= 0 || space_penalty == 0.0) {
			vector<DiagGraph::Edge>::const_iterator edge = diags.get_edge(d_idx, d.j);
			if (edge != diags.edges.end() && edge->prefix_score > e.prefix_score + gap_score + d.score)
				return 0;
			/*if (d.prefix_score > e.prefix_score + gap_score + d.score)
				return 0;*/
			Link link;
			if (get_link(e, d, query, subject, link, link_padding) > 0) {
				diff1 = e.score - link.score1;
				const int prefix_e = diags.prefix_score(e_idx, link.subject_pos1, path_max, path_min);
				prefix_score = prefix_e - diff1 + gap_score + link.score2;
				vector<DiagGraph::Edge>::const_iterator edge = diags.get_edge(d_idx, link.subject_pos2);
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
			vector<DiagGraph::Edge>::const_iterator edge = diags.get_edge(d_idx, d.j);
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
			diags.add_edge(DiagGraph::Edge(prefix_score, path_max, link_j, d_idx, e_idx, prefix_score == path_max ? prefix_score : path_min, prefix_score_begin));
			if (log)
				cout << "Link n=" << e_idx << " d=" << e.diag() << " i_end=" << e.query_end() << " max_i=" << max_i << " shift=" << shift << " space=" << space << " prefix_score=" << prefix_score << " link_score=" << link_score << " path_min="<<path_min<<endl;
		}
		return prefix_score;
	}

	template<typename _it>
	void Aligner::forward_pass(_it begin, _it end, bool init, double space_penalty)
	{
		window.clear();

		for (_it it = begin; it != end; ++it) {

			unsigned node = (unsigned)(*it);
			if(init) diags.init(node);
			DiagonalNode& d = diags[node];
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
				DiagonalNode&e = diags[j->second];
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
				DiagonalNode &e = diags[j->second];
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

	int Aligner::run(list<Hsp> &hsps, list<ApproxHsp> &ts, double space_penalty, int cutoff, int max_shift)
	{
		if (config.chaining_maxnodes > 0) {
			std::sort(diags.nodes.begin(), diags.nodes.end(), DiagonalSegment::cmp_score);
			if (diags.nodes.size() > config.chaining_maxnodes)
				diags.nodes.erase(diags.nodes.begin() + config.chaining_maxnodes, diags.nodes.end());
		}
		if (config.chaining_len_cap > 0.0 && diags.nodes.size() > config.chaining_min_nodes) {
			std::sort(diags.nodes.begin(), diags.nodes.end(), DiagonalSegment::cmp_score);
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

		forward_pass(IndexIterator(0llu), IndexIterator(diags.nodes.size()), true, space_penalty);
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

	int Aligner::run(list<Hsp> &hsps, list<ApproxHsp> &ts, vector<DiagonalSegment>::const_iterator begin, vector<DiagonalSegment>::const_iterator end, int band)
	{
		if (log)
			cout << "***** Seed hit run " << begin->diag() << '\t' << (end - 1)->diag() << '\t' << (end - 1)->diag() - begin->diag() << endl;
		diags.init();
		diags.load(begin, end);
		return run(hsps, ts, SPACE_PENALTY, 19, band);
	}

	Aligner::Aligner(const Sequence &query, const Sequence &subject, bool log, unsigned frame) :
		query(query),
		subject(subject),
		//query_bc(query_bc),
		log(log),
		frame(frame)
	{
	}

thread_local DiagGraph Aligner::diags;
thread_local map<int, unsigned> Aligner::window;

}

namespace Chaining {

static Score merge_score(const ApproxHsp& h1, const ApproxHsp& h2) {
	static const double GAP_PENALTY = 0.5;
	const Loc gq = h2.query_range.begin_ - h1.query_range.end_, gt = h2.subject_range.begin_ - h1.subject_range.end_;
	if (gq < 0 || gt < 0)
		return 0;
	const Score s = h1.score + h2.score;
	if (gq > gt) {
		return Score(s - gq * GAP_PENALTY - gt * SPACE_PENALTY);
	}
	else
		return Score(s - gt * GAP_PENALTY - gq * SPACE_PENALTY);
}

static ApproxHsp merge(const ApproxHsp& h1, const ApproxHsp& h2) {
	ApproxHsp h(h1.frame);
	h.d_max = max(h1.d_max, h2.d_max);
	h.d_min = min(h1.d_min, h2.d_min);
	h.query_range = { h1.query_range.begin_, h2.query_range.end_ };
	h.query_source_range = h.query_range;
	h.subject_range = { h1.subject_range.begin_, h2.subject_range.end_ };
	h.score = merge_score(h1, h2);
	h.evalue = 0;
	if (h1.max_diag.score > h2.max_diag.score) {
		h.max_diag = h1.max_diag;
		h.max_diag.d_max_right = max(h.max_diag.d_max_right, h2.d_max);
		h.max_diag.d_min_right = min(h.max_diag.d_min_right, h2.d_min);
	}
	else {
		h.max_diag = h2.max_diag;
		h.max_diag.d_max_left = max(h.max_diag.d_max_left, h1.d_max);
		h.max_diag.d_min_left = min(h.max_diag.d_min_left, h1.d_min);
	}
	return h;
}

static void merge_hsps(list<ApproxHsp>& hsps) {
	auto it = hsps.begin();
	while (it != hsps.end()) {
		auto it2 = it;
		++it2;
		while (it2 != hsps.end()) {
			if (merge_score(*it, *it2) > max(it->score, it2->score)) {
				*it = merge(*it, *it2);
				it2 = hsps.erase(it2);
			}
			else if (merge_score(*it2, *it) > max(it->score, it2->score)) {
				*it = merge(*it2, *it);
				it2 = hsps.erase(it2);
			}
			else
				++it2;
		}
		++it;
	}
}

std::pair<int, list<ApproxHsp>> run(Sequence query, Sequence subject, vector<DiagonalSegment>::const_iterator begin, vector<DiagonalSegment>::const_iterator end, bool log, unsigned frame)
{
	const int band = config.chaining_maxgap;
	if (end - begin == 1) {
		const Loc d = begin->diag();
		const Anchor anchor(*begin, d, d, d, d, begin->score);
		return { begin->score, { { d, d, begin->score, (int)frame, begin->query_range(), begin->subject_range(), anchor}} };
	}
	Chaining::Aligner ga(query, subject, log, frame);
	list<Hsp> hsps;
	list<ApproxHsp> ts;
	int score = ga.run(hsps, ts, begin, end, band);
	if (!config.no_chaining_merge_hsps)
		merge_hsps(ts);
	return std::make_pair(score, std::move(ts));
}

list<Hsp> run(Sequence query, const std::vector<DpTarget>& targets) {
	list<Hsp> out;
	return out;
}

}