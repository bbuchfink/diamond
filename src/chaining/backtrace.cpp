#include <iostream>
#include "aligner.h"
#include "../basic/config.h"
#include "../util/geo/geo.h"

using std::vector;
using std::cout;
using std::endl;
using std::list;
using std::min;
using std::max;
using std::numeric_limits;

namespace Chaining {

static bool disjoint(list<ApproxHsp>::const_iterator begin, list<ApproxHsp>::const_iterator end, const ApproxHsp &t, int cutoff)
{
	for (; begin != end; ++begin) {
		const double ot = t.subject_range.overlap_factor(begin->subject_range),
			oq = t.query_range.overlap_factor(begin->query_range);
		if ((1.0 - min(ot, oq)) * t.score / begin->score >= config.chaining_stacked_hsp_ratio)
			continue;
		if ((1.0 - max(ot, oq)) * t.score < cutoff)
			return false;
		//if (begin->partial_score(t) < cutoff || !begin->collinear(t))
		//if (!begin->disjoint(t) || !begin->collinear(t))
		//if (!begin->rel_disjoint(t))
		//	return false;
	}
	return true;
}

static bool disjoint(list<ApproxHsp>::const_iterator begin, list<ApproxHsp>::const_iterator end, const DiagonalSegment &d, int cutoff)
{
	for (; begin != end; ++begin) {
		const double ot = d.subject_range().overlap_factor(begin->subject_range),
			oq = d.query_range().overlap_factor(begin->query_range);
		if ((1.0 - min(ot, oq)) * d.score / begin->score >= config.chaining_stacked_hsp_ratio)
			continue;
		if ((1.0 - max(ot, oq)) * d.score < cutoff)
			return false;
		//if (begin->partial_score(d) < cutoff || !begin->collinear(d))
		//if (!begin->disjoint(d) || !begin->collinear(d))
		//if (!begin->rel_disjoint(d))
		//return false;
	}
	return true;
}

bool Aligner::backtrace_old(size_t node, int j_end, Hsp* out, ApproxHsp& t, int score_max, int score_min, int max_shift, unsigned& next) const
{
	
		const DiagonalNode& d = diags[node];
		vector<DiagGraph::Edge>::const_iterator f = diags.get_edge(node, j_end);
		bool at_end = f >= diags.edges.end();
		const int prefix_score = at_end ? d.score : f->prefix_score;
		if (prefix_score > score_max)
			return false;

		int j;
		score_min = std::min(score_min, at_end ? 0 : f->prefix_score_begin);

		//if (f != diags.edges.end() && (!stop_at_min || f->path_min == diags[f->node_out].path_min)) {
		if (!at_end) {
			const DiagonalNode& e = diags[f->node_out];
			const int shift = d.diag() - e.diag();
			j = f->j;

			if (abs(shift) <= max_shift) {
				const bool bt = backtrace_old(f->node_out, shift > 0 ? j : j + shift, out, t, score_max, score_min, max_shift, next);
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
			const DiagonalNode& e = diags[f->node_out];
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
		if (d.score > t.max_diag.score) {
			t.max_diag = d;
			t.max_diag.prefix_score = prefix_score;
			t.max_diag.d_max_left = max(max(t.max_diag.d_max_right, t.max_diag.d_max_left), dd);
			t.max_diag.d_min_left = min(min(t.max_diag.d_min_right, t.max_diag.d_min_left), dd);
			t.max_diag.d_max_right = dd;
			t.max_diag.d_min_right = dd;
		}
		else {
			t.max_diag.d_max_right = std::max(t.max_diag.d_max_right, dd);
			t.max_diag.d_min_right = std::min(t.max_diag.d_min_right, dd);
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
		return true;

}



void Aligner::backtrace(const size_t node, const int j_end, Hsp* out, ApproxHsp& t, const int score_max, const int score_min, const int max_shift, unsigned& next) const
{
	struct Node {
		size_t node;
		int score_min, j_end;
	};
	vector<Node> nodes;
	nodes.push_back({ node, score_min, j_end });
	bool ret = false, ret_value;

	while (!nodes.empty()) {
		Node& c = nodes.back();
		const DiagonalNode& d = diags[c.node];
		vector<DiagGraph::Edge>::const_iterator f = diags.get_edge(c.node, c.j_end);
		bool at_end = f >= diags.edges.end();
		const int prefix_score = at_end ? d.score : f->prefix_score;
		if (!ret && prefix_score > score_max) {
			ret = true;
			ret_value = false;
			nodes.pop_back();
			continue;
		}

		c.score_min = std::min(c.score_min, at_end ? 0 : f->prefix_score_begin);

		//if (f != diags.edges.end() && (!stop_at_min || f->path_min == diags[f->node_out].path_min)) {
		if (!at_end) {
			const DiagonalNode& e = diags[f->node_out];
			const int shift = d.diag() - e.diag();
			const int j = f->j;

			if (abs(shift) <= max_shift) {
				if (!ret) {
					nodes.push_back({ f->node_out, c.score_min, shift > 0 ? j : j + shift });
					continue;
				}
				assert(ret == true);
				if (!ret_value) {
					if (f->prefix_score_begin > c.score_min) {
						ret_value = false;
						nodes.pop_back();
						continue;
					}
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
				out->score = score_max - c.score_min;
			}
			t.query_range.begin_ = d.i;
			t.subject_range.begin_ = d.j;
			t.score = score_max - c.score_min;
		}
		else {
			const DiagonalNode& e = diags[f->node_out];
			const int shift = d.diag() - e.diag();
			if (out) {
				if (shift > 0) {
					out->transcript.push_back(op_insertion, (unsigned)shift);
					out->length += shift;
				}
				else if (shift < 0) {
					const int j = f->j;
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
			int j = at_end ? d.j : f->j;
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
		ret = true;
		ret_value = true;
		nodes.pop_back();
	}
}

void Aligner::backtrace(size_t top_node, Hsp* out, ApproxHsp& t, int max_shift, unsigned& next, int max_j) const
{
	ApproxHsp traits(frame);
	if (top_node != DiagGraph::end) {
		const DiagonalNode& d = diags[top_node];
		if (out) {
			out->transcript.clear();
			out->query_range.end_ = d.query_end();
			out->subject_range.end_ = d.subject_end();
		}
		traits.subject_range.end_ = d.subject_end();
		traits.query_range.end_ = d.query_end();
		int score_min = d.prefix_score;
		backtrace_old(top_node, std::min(d.subject_end(), max_j), out, traits, d.prefix_score, score_min, max_shift, next);
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

int Aligner::backtrace(size_t top_node, list<Hsp>& hsps, list<ApproxHsp>& ts, list<ApproxHsp>::iterator& t_begin, int cutoff, int max_shift) const
{
	unsigned next;
	int max_score = 0, max_j = (int)subject.length();
	do {
		Hsp* hsp = log ? new Hsp(true) : 0;
		ApproxHsp t(frame);
		next = std::numeric_limits<unsigned>::max();
		backtrace(top_node, hsp, t, max_shift, next, max_j);
		//assert(t.score > 0);
		//Geo::assert_diag_bounds(t.d_max, query.length(), subject.length());
		//Geo::assert_diag_bounds(t.d_min, query.length(), subject.length());
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

int Aligner::backtrace(list<Hsp>& hsps, list<ApproxHsp>& ts, int cutoff, int max_shift) const
{
	vector<DiagonalNode*> top_nodes;
	for (size_t i = 0; i < diags.nodes.size(); ++i) {
		DiagonalNode& d = diags.nodes[i];
		//cout << "node=" << i << " prefix_score=" << d.prefix_score << " path_max=" << d.path_max << " rel_score=" << d.rel_score() << " cutoff=" << cutoff << endl;
		//if (d.prefix_score >= cutoff && (d.prefix_score == d.path_max || d.prefix_score - d.path_min >= cutoff))
		if (d.rel_score() >= cutoff)
			top_nodes.push_back(&d);
	}
	std::sort(top_nodes.begin(), top_nodes.end(), DiagonalNode::cmp_rel_score);
	int max_score = 0;
	list<ApproxHsp>::iterator t_begin = ts.end();

	for (vector<DiagonalNode*>::const_iterator i = top_nodes.begin(); i < top_nodes.end(); ++i) {
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

}