#pragma once
#include "../util/geo/diagonal_segment.h"

struct DiagonalNode : public DiagonalSegment
{
	enum { estimate, finished };
	DiagonalNode() :
		DiagonalSegment(),
		link_idx(-1),
		prefix_score(0),
		path_max(0),
		path_min(0)
	{}
	DiagonalNode(int query_pos, int subject_pos, int len, int score, int link_idx = -1) :
		DiagonalSegment(query_pos, subject_pos, len, score),
		link_idx(link_idx),
		prefix_score(score),
		path_max(score),
		path_min(score)
	{}
	DiagonalNode(const DiagonalSegment& d) :
		DiagonalSegment(d),
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
	static bool cmp_prefix_score(const DiagonalNode* x, const DiagonalNode* y)
	{
		return x->prefix_score > y->prefix_score;
	}
	static bool cmp_rel_score(const DiagonalNode* x, const DiagonalNode* y)
	{
		return x->rel_score() > y->rel_score();
	}
	int link_idx, prefix_score, path_max, path_min;
};

struct DiagGraph
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

	void load(std::vector<DiagonalSegment>::const_iterator begin, std::vector<DiagonalSegment>::const_iterator end);
	void sort();
	void prune();
	void clear_edges();

	std::vector<Edge>::iterator add_edge(const Edge& edge)
	{
		for (std::vector<DiagonalNode>::iterator j = nodes.begin() + edge.node_in + 1; j < nodes.end(); ++j)
			if (j->link_idx == -1)
				break;
			else
				++j->link_idx;
		assert(nodes[edge.node_in].link_idx >= 0 && nodes[edge.node_in].link_idx <= (int)edges.size());
		DiagonalNode& d = nodes[edge.node_in];
		if (edge.prefix_score > d.prefix_score) {
			d.prefix_score = edge.prefix_score;
			d.path_max = edge.path_max;
			d.path_min = edge.path_min;
		}
		return edges.insert(edges.begin() + d.link_idx++, edge);
	}

	std::vector<Edge>::const_iterator get_edge(size_t node, int j) const
	{
		const DiagonalNode& d = nodes[node];
		if (d.score == 0)
			return edges.begin() + d.link_idx - 1;
		if (edges.empty())
			return edges.end();
		int max_score = d.score;
		ptrdiff_t max_i = -1;
		for (ptrdiff_t i = d.link_idx - 1; i >= 0 && edges[i].node_in == node; --i)
			if (edges[i].j < j && edges[i].prefix_score > max_score) {
				max_i = i;
				max_score = edges[i].prefix_score;
			}

		/*for (vector<Edge>::const_iterator i = edges.begin() + d.link_idx - 1; i >= edges.begin() && i->node_in == node; --i)
			if (i->j < j && i->prefix_score > max_score) {
				max_edge = i;
				max_score = i->prefix_score;
			}*/
		return max_i >= 0 ? edges.begin() + max_i : edges.end();
	}

	int prefix_score(size_t node, int j, int& path_max, int& path_min) const
	{
		const std::vector<Edge>::const_iterator i = get_edge(node, j);
		path_max = i == edges.end() ? nodes[node].score : std::max(nodes[node].score, i->path_max);
		path_min = i == edges.end() ? nodes[node].score : i->path_min;
		return i == edges.end() ? nodes[node].score : std::max(nodes[node].score, i->prefix_score);
	}

	DiagonalNode& operator[](size_t k)
	{
		return nodes[k];
	}

	const DiagonalNode& operator[](size_t k) const
	{
		return nodes[k];
	}

	void print(Sequence query, Sequence subject) const;
	size_t top_node() const;

	std::vector<DiagonalNode> nodes;
	std::vector<Edge> edges;
};

void smith_waterman(Sequence q, Sequence s, const DiagGraph& diags);