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

	enum { link_padding = 5, reverse_link_min_overhang = 10 };

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
		
	int follow_path(unsigned level, unsigned node, int score, int subject_pos)
	{
		static const int min_nw_space = 7;
		static const float space_penalty = -0.5;
		Diagonal_node* d = &diags[node];
		if (log) indent(cout, level) << "Node " << node << " Score=" << d->prefix_score() << " j=" << subject_pos << endl;
		int max_score = d->score;

		for (unsigned k = 0; k < Diagonal_node::n_path; ++k) {
			Diagonal_node::Edge *f = &d->edges[k];
			if (f->prefix_score == 0)
				break;
			if (f->j > subject_pos)
				continue;
			if (f->state == Diagonal_node::finished) {
				if (log)
					indent(cout, level) << "Visited node final_score=" << f->prefix_score << endl;
				return f->prefix_score;
			}
			unsigned next = f->node;
			const Diagonal_node *e = &diags[next];
			const int shift = d->diag() - e->diag();
			const int estimate_score = f->prefix_score;
			if (!f->exact) {
				const int space = shift > 0 ? d->j - e->subject_last() : d->i - e->query_last();
				if (log)
					indent(cout, level) << "Node=" << node << " Link n=" << f->node << " space=" << space << " shift=" << shift << endl;
				if (space >= min_nw_space && abs(shift) > 1) {
					if (log) {
						const sequence q1 = sequence(query, e->query_last() + 1, d->i - 1),
							s1 = sequence(subject, e->subject_last() + 1, d->j - 1);
						indent(cout, level) << q1 << endl;
						indent(cout, level) << s1 << endl;
					}
					const int nw = needleman_wunsch(query, subject, e->query_last() + 1, d->i, e->subject_last() + 1, d->j, node, k, diags, log);
					d = &diags[node];
					f = &d->edges[k];
					f->prefix_score = nw + d->score;
				}
				else {
					Link l;
					int link_score = get_link(*e, *d, query, subject, l, link_padding);
					f->prefix_score = l.score1 - e->score - config.gap_open - abs(shift)*config.gap_extend + l.score2;
					f->j = l.subject_pos2;
					if (log)
						indent(cout, level) << "Direct link scores=" << l.score1 << "," << l.score2 << endl;
				}
			}
			else
				f->prefix_score = -f->diff1 - config.gap_open - abs(shift)*config.gap_extend + d->score - f->diff2;
			if (log)
				indent(cout, level) << "Node=" << node << " Link n=" << next << " estimate_score=" << estimate_score << " prefix_score=" << f->prefix_score << " j=" << f->j << endl;
			const int fp_score = follow_path(level + 1, next, 0, f->j + std::min(shift, 0) - 1);
			d = &diags[node];
			f = &d->edges[k];
			f->prefix_score += fp_score;
			f->state = Diagonal_node::finished;
			if (log)
				indent(cout, level) << "Node=" << node << " Link n=" << next << " final prefix_score=" << f->prefix_score << endl;
			max_score = std::max(max_score, f->prefix_score);
		}

		top_node = std::min(top_node, Node_ref(node, max_score));
		d->edges.sort();
		if(log) indent(cout, level) << "Node " << node << " final_score=" << max_score << endl;
		return max_score;
	}

	int get_approximate_link(Diagonal_node &d, const Diagonal_node &e, int k)
	{
		const int shift = d.diag() - e.diag();
		int gap_score = -config.gap_open - abs(shift)*config.gap_extend;
		const int space = shift > 0 ? d.j - e.subject_last() : d.i - e.query_last();
		int prefix_score = 0, link_score = 0, link_j, diff1 = 0, diff2 = 0;
		bool exact;
		if (space <= 0) {
			Link link;
			if (get_link(e, d, query, subject, link, link_padding) > 0) {
				diff1 = e.score - link.score1;
				diff2 = d.score - link.score2;
				prefix_score = e.prefix_score(link.subject_pos1) - diff1 + gap_score + link.score2;
				link_score = link.score1 + link.score2 + gap_score;
				exact = true;
				link_j = link.subject_pos2;
			}
		}
		else {
			prefix_score = e.prefix_score(e.subject_last()) + gap_score - int(config.raw_space_penalty*std::max(space - 1, 0)) + d.score;
			link_score = e.score + d.score + gap_score;
			exact = false;
			link_j = d.j;
		}

		if (log)
			cout << "Link n=" << k << " shift=" << shift << " space=" << space << " prefix_score=" << prefix_score << " link_score=" << link_score << endl;
		if (prefix_score > 0)
			d.edges.add(Diagonal_node::Edge(prefix_score, link_j, k, exact, Diagonal_node::estimate, diff1, diff2));
		return prefix_score;
	}

	void forward_pass()
	{
		static const int max_dist = 60;
		static const float space_penalty = -0.5;
		
		for (unsigned node = 0; node < diags.size(); ++node) {
			Diagonal_node& d = diags[node];
			if (log) cout << "Node " << node << " Score=" << d.score << endl;
			for (int k = node - 1; k >= 0; --k) {
				const Diagonal_node &e = diags[k];
				if (std::max(d.j - e.subject_last(), 0) < max_dist) {
					if (abs(d.diag() - e.diag()) > max_dist)
						continue;
					get_approximate_link(d, e, k);
				}
				else
					break;
			}
			if (log)
				cout << "Final score=" << d.prefix_score() << endl;
			top_nodes.push_back(Node_ref(node, d.prefix_score()));

			for (int k = 0; k < Diagonal_node::n_path; ++k) {
				Diagonal_node::Edge &f = d.edges[k];
				if (f.prefix_score == 0)
					break;
				Diagonal_node &e = diags[f.node];
				//if (e.subject_end() - (d.subject_end() - std::min(e.diag()-d.diag(), 0)) >= reverse_link_min_overhang) {
				if (e.j + std::min(e.diag() - d.diag(), 0) - d.j >= reverse_link_min_overhang) {
					if (log)
						cout << "Computing reverse link node=" << f.node << endl;
					const int s = e.prefix_score(f.j + std::min(d.diag() - e.diag(), 0));
					const int s2 = get_approximate_link(e, d, node);
					f.prefix_score += s2 - s;
					if (log)
						cout << "New prefix score=" << f.prefix_score << endl;
				}
			}
			if (log)
				cout << endl;
		}
	}

	void backtrace(unsigned node, int j_end, Hsp_data &out)
	{
		const Diagonal_node &d = diags[node];
		int j, edge = d.get_edge(j_end);
		if (edge >= 0) {
			const Diagonal_node::Edge &f = d.edges[edge];
			const Diagonal_node &e = diags[f.node];
			const int shift = d.diag() - e.diag();
			j = f.j;

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
			j = d.j;
		}

		const int d2 = d.diag();
		if (log) cout << "Backtrace node=" << node << " i=" << d2 + j << "-" << d2 + j_end << " j=" << j << "-" << j_end << endl;
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
		diag_scores(TLS::get(diag_scores_ptr)),
		diags(TLS::get(diags_ptr)),
		top_nodes(TLS::get(top_nodes_ptr))
	{
		static const double score_range = 0.97;

		config.min_diag_raw_score = 15;
		diags.clear();
		top_nodes.clear();
		diag_scores.scan_diags(sh[0], query, subject, qp, log, diags);
		std::sort(diags.begin(), diags.end(), Diagonal_segment::cmp_heuristic);
		if (log)
			for (int k = 0; k < (int)diags.size(); ++k) {
				const Diagonal_segment &d = diags[k];
				cout << "Diag n=" << k << " i=" << d.i << " j=" << d.j << " score=" << d.score << " len=" << d.len << endl;
				cout << sequence(query, d.i, d.query_last()) << endl;
				cout << sequence(subject, d.j, d.subject_last()) << endl;
			}
		if (log) cout << endl << endl;

		forward_pass();

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
			cout << buf.get_begin() << endl << "Smith-Waterman:" << endl;
			smith_waterman(query, subject, diags);
			cout << endl << endl;
		}
	}

	static TLS_PTR Diag_scores *diag_scores_ptr;
	static TLS_PTR vector<Diagonal_node> *diags_ptr;
	static TLS_PTR vector<Node_ref> *top_nodes_ptr;
	const sequence query, subject;
	const Long_score_profile &qp;
	const bool log;
	
	Diag_scores &diag_scores;
	vector<Diagonal_node> &diags;
	vector<Node_ref> &top_nodes;
	Node_ref top_node;

};

TLS_PTR Diag_scores *Greedy_aligner2::diag_scores_ptr;
TLS_PTR vector<Diagonal_node> *Greedy_aligner2::diags_ptr;
TLS_PTR vector<Greedy_aligner2::Node_ref> *Greedy_aligner2::top_nodes_ptr;

void greedy_align2(sequence query, const Long_score_profile &qp, sequence subject, const vector<Diagonal_segment> &sh, bool log, Hsp_data &out)
{
	Greedy_aligner2(query, qp, subject, sh, log, out);
}