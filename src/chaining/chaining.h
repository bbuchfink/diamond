/****
DIAMOND protein aligner
Copyright (C) 2013-2019 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef CHAINING_H_
#define CHAINING_H_

#include <utility>
#include <list>
#include "../dp/hsp_traits.h"
#include "../dp/comp_based_stats.h"

std::pair<int, std::list<Hsp_traits>> greedy_align(sequence query, sequence subject, std::vector<Diagonal_segment>::const_iterator begin, std::vector<Diagonal_segment>::const_iterator end, bool log, unsigned frame);

struct Diagonal_node : public Diagonal_segment
{
	enum { estimate, finished };
	Diagonal_node() :
		Diagonal_segment(),
		link_idx(-1),
		prefix_score(0),
		path_max(0),
		path_min(0)
	{}
	Diagonal_node(int query_pos, int subject_pos, int len, int score, int link_idx = -1) :
		Diagonal_segment(query_pos, subject_pos, len, score),
		link_idx(link_idx),
		prefix_score(score),
		path_max(score),
		path_min(score)
	{}
	Diagonal_node(const Diagonal_segment &d) :
		Diagonal_segment(d),
		link_idx(-1),
		prefix_score(d.score),
		path_max(d.score),
		path_min(d.score)
	{}
	void deactivate()
	{
		link_idx = 0;
	}
	void reset()
	{
		link_idx = -1;
		prefix_score = score;
		path_max = score;
		path_min = score;
	}
	bool is_maximum() const
	{
		return path_max == prefix_score;
	}
	int rel_score() const
	{
		return prefix_score == path_max ? prefix_score : prefix_score - path_min;
	}
	static bool cmp_prefix_score(const Diagonal_node *x, const Diagonal_node *y)
	{
		return x->prefix_score > y->prefix_score;
	}
	static bool cmp_rel_score(const Diagonal_node *x, const Diagonal_node *y)
	{
		return x->rel_score() > y->rel_score();
	}
	int link_idx, prefix_score, path_max, path_min;
};

struct Diag_graph
{

	enum { end = size_t(-1) };

	struct Edge
	{
		Edge() :
			prefix_score(0),
			node_in()
		{
		}
		Edge(int prefix_score, int path_max, int j, unsigned node_in, unsigned node_out, int path_min, int prefix_score_begin) :
			prefix_score(prefix_score),
			path_max(path_max),
			j(j),
			path_min(path_min),
			prefix_score_begin(prefix_score_begin),
			node_in(node_in),
			node_out(node_out)
		{
		}
		/*operator int() const
		{
			return prefix_score;
		}
		bool operator<(const Edge &x) const
		{
			return prefix_score > x.prefix_score;
		}*/
		int prefix_score, path_max, j, path_min, prefix_score_begin;
		unsigned node_in, node_out;
	};

	void init()
	{
		nodes.clear();
		edges.clear();
	}

	void init(unsigned node)
	{
		if (edges.size() >= (size_t)std::numeric_limits<int>::max())
			throw std::runtime_error("Too many edges.");
		nodes[node].link_idx = (int)edges.size();
	}

	void load(vector<Diagonal_segment>::const_iterator begin, vector<Diagonal_segment>::const_iterator end);
	void sort();
	void prune();
	void clear_edges();

	vector<Edge>::iterator add_edge(const Edge &edge)
	{
		for (vector<Diagonal_node>::iterator j = nodes.begin() + edge.node_in + 1; j < nodes.end(); ++j)
			if (j->link_idx == -1)
				break;
			else
				++j->link_idx;
		assert(nodes[edge.node_in].link_idx >= 0 && nodes[edge.node_in].link_idx <= (int)edges.size());
		Diagonal_node &d = nodes[edge.node_in];
		if (edge.prefix_score > d.prefix_score) {
			d.prefix_score = edge.prefix_score;
			d.path_max = edge.path_max;
			d.path_min = edge.path_min;
		}
		return edges.insert(edges.begin() + d.link_idx++, edge);
	}

	vector<Edge>::const_iterator get_edge(size_t node, int j) const
	{
		const Diagonal_node &d = nodes[node];
		if (d.score == 0)
			return edges.begin() + d.link_idx - 1;
		if (edges.empty())
			return edges.end();
		int max_score = d.score;
		vector<Edge>::const_iterator max_edge = edges.end();
		for (vector<Edge>::const_iterator i = edges.begin() + d.link_idx - 1; i >= edges.begin() && i->node_in == node; --i)
			if (i->j < j && i->prefix_score > max_score) {
				max_edge = i;
				max_score = i->prefix_score;
			}
		return max_edge;
	}

	int prefix_score(size_t node, int j, int &path_max, int &path_min) const
	{
		const vector<Edge>::const_iterator i = get_edge(node, j);
		path_max = i == edges.end() ? nodes[node].score : std::max(nodes[node].score, i->path_max);
		path_min = i == edges.end() ? nodes[node].score : i->path_min;
		return i == edges.end() ? nodes[node].score : std::max(nodes[node].score, i->prefix_score);
	}

	Diagonal_node& operator[](size_t k)
	{
		return nodes[k];
	}

	const Diagonal_node& operator[](size_t k) const
	{
		return nodes[k];
	}

	void print(sequence query, sequence subject) const;
	size_t top_node() const;

	vector<Diagonal_node> nodes;
	vector<Edge> edges;
};

void smith_waterman(sequence q, sequence s, const Diag_graph &diags);

#endif