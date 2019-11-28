#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <limits.h>
#include <stdlib.h>
#include <iomanip>
#include <unordered_map>
#include "../io/text_input_file.h"
#include "../../basic/config.h"
#include "../string/tokenizer.h"
#include "../log_stream.h"
#include "../../lib/MemoryPool/MemoryPool.h"
#include "edge_vec.h"

using std::string;
using std::list;
using std::endl;
using std::map;
using std::vector;
using std::ostream;
using std::pair;
using std::cout;

namespace Util { namespace Algo { namespace UPGMA_MC {

struct Edge {
	Edge(int n1, int n2, int count, double s):
		n1(n1),
		n2(n2),
		deleted(0),
		count(count),
		s(s),
		l(s),
		u(s)
	{}
	void set_bounds(double lambda, double max_dist, double max_edges) {
		l = (s + lambda * (max_edges - count)) / max_edges;
		u = (s + max_dist * (max_edges - count)) / max_edges;
	}
	bool operator<=(const Edge &e) const {
		return u <= e.l;
	}
	int target(int me) const {
		if (n1 == me)
			return n2;
		else if (n2 == me)
			return n1;
		else
			throw std::runtime_error("Edge::target");
	}
	friend ostream& operator<<(ostream &s, const Edge &e) {
		s << e.n1 << '\t' << e.n2 << '\t' << e.l << '\t' << e.u;
		return s;
	}
	int n1, n2, deleted, count;
	double s, l, u;
};

typedef list<Edge, MemoryPool<Edge>> EdgeList;
typedef EdgeList::iterator EdgePtr;

struct CmpEdge {
	bool operator()(const EdgePtr& e, const EdgePtr &f) const {
		return e->l > f->l || (e->l == f->l && (e->n1 > f->n1 || (e->n1 == f->n1 && e->n2 > f->n2)));
	}
};

typedef std::priority_queue<EdgePtr, vector<EdgePtr>, CmpEdge> Queue;

void erase(EdgePtr &e, EdgeList &edges) {
	++e->deleted;
	if (e->deleted == 3)
		edges.erase(e);
}

struct Node {
	Node(int idx, int size, int parent):
		idx(idx),
		size(size),
		parent(parent)
	{}
	void sort_neighbors() {
		std::sort(neighbors.begin(), neighbors.end(), CmpNeighbor{ parent });
	}
	void set_parent(int parent, EdgeList &edges) {
		this->parent = parent;
		for (EdgePtr e : neighbors)
			erase(e, edges);
		neighbors.clear();
		neighbors.shrink_to_fit();
	}
	bool root() const {
		return parent == idx;
	}
	struct CmpNeighbor {
		int me;
		bool operator()(const EdgePtr &e, const EdgePtr &f) const { 
			return e->target(me) < f->target(me);
		}
	};
	int idx, size, parent;
	vector<EdgePtr> neighbors;
};

bool valid(const EdgePtr &e, const vector<Node> &nodes) {
	return nodes[e->n1].root() && nodes[e->n2].root();
}

void merge_nodes(int n1,
	int n2,
	vector<Node> &nodes,
	EdgeList &edges,
	Queue &queue,
	double max_dist,
	double lambda) {	
	const int union_idx = (int)nodes.size();
	nodes.emplace_back(union_idx, nodes[n1].size + nodes[n2].size, union_idx);
	Node &node1 = nodes[n1], &node2 = nodes[n2], &union_node = nodes.back();
	
	vector<EdgePtr>::iterator i = node1.neighbors.begin(), j = node2.neighbors.begin();
	while (i < node1.neighbors.end() || j < node2.neighbors.end()) {
		int it = i < node1.neighbors.end() ? (*i)->target(n1) : INT_MAX, jt = j < node2.neighbors.end() ? (*j)->target(n2) : INT_MAX;
		double s;
		int edge_count;
		if (it == jt) {
			s = (*i)->s + (*j)->s;
			edge_count = (*i)->count + (*j)->count;
			++i;
			++j;
		}
		else if (it < jt) {
			s = (*i)->s;
			edge_count = (*i)->count;
			++i;
		}
		else {
			s = (*j)->s;
			edge_count = (*j)->count;
			it = jt;
			++j;
		}
		if (nodes[it].parent != it || it == n1 || it == n2)
			continue;
		const double max_edges = (double)union_node.size*(double)nodes[it].size;		
		const EdgePtr e = edges.emplace(edges.end(), it, union_idx, edge_count, s);
		e->set_bounds(lambda, max_dist, max_edges);
		queue.push(e);
		union_node.neighbors.push_back(e);
		nodes[it].neighbors.push_back(e);
	}
	node1.set_parent(union_idx, edges);
	node2.set_parent(union_idx, edges);
}

int parent(int idx, const vector<Node> &nodes) {
	while (nodes[idx].parent != idx)
		idx = nodes[idx].parent;
	return idx;
}

struct PairHash {
	size_t operator()(const pair<int, int> &x) const {
		return std::hash<int>()(x.first) ^ std::hash<int>()(x.second);
	}
};

double load_edges(EdgeVec::const_iterator &begin, const EdgeVec::const_iterator &end, EdgeList &edges, vector<Node> &nodes, Queue &queue, double lambda, double max_dist) {
	message_stream << "Clearing neighborhoods..." << endl;
	for (Node &node : nodes) {
		node.neighbors.clear();
		node.neighbors.shrink_to_fit();
	}

	message_stream << "Clearing old edges..." << endl;
	for (EdgeList::iterator i = edges.begin(); i != edges.end();)
		if (!valid(i, nodes))
			i = edges.erase(i);
		else {
			i->deleted = 0;
			++i;
		}

	if (edges.size() >= config.upgma_edge_limit)
		throw std::runtime_error("Edge limit");

	message_stream << "Building edge hash map..." << endl;
	std::unordered_map<pair<int, int>, EdgePtr, PairHash> edge_map;
	edge_map.reserve(edges.size());
	for (EdgePtr i = edges.begin(); i != edges.end(); ++i)
		edge_map[{i->n1, i->n2}] = i;

	double evalue = lambda;
	message_stream << "Reading edges..." << endl;
	while (edges.size() < config.upgma_edge_limit && begin < end) {
		const int query_idx = begin->n1, target_idx = begin->n2;
		evalue = begin->d;

		int i = parent(query_idx, nodes), j = parent(target_idx, nodes);
		if (i == query_idx && j == target_idx) {
			if (i >= j) std::swap(i, j);
			edges.emplace_back(i, j, 1, evalue);
		}
		else {
			if (i >= j) std::swap(i, j);
			const auto e = edge_map.find({ i,j });
			if (e == edge_map.end()) {
				const EdgePtr f = edges.emplace(edges.end(), i, j, 1, evalue);
				edge_map[{i, j}] = f;
			}
			else {
				++e->second->count;
				e->second->s += evalue;
			}
		}
	
		++begin;
		/*if (edges.size() % 10000 == 0 && !edges.empty())
			message_stream << "#Edges: " << edges.size() << endl;*/
	}

	lambda = (begin == end) ? max_dist : evalue;

	message_stream << "Recomputing bounds, building edge vector and neighborhood..." << endl;
	vector<EdgePtr> edge_vec;
	edge_vec.reserve(edges.size());
	for (EdgeList::iterator i = edges.begin(); i != edges.end(); ++i) {
		i->set_bounds(lambda, max_dist, (double)nodes[i->n1].size * (double)nodes[i->n2].size);
		edge_vec.push_back(i);
		nodes[i->n1].neighbors.push_back(i);
		nodes[i->n2].neighbors.push_back(i);
	}

	message_stream << "Sorting neighborhoods..." << endl;
	for (Node &node : nodes)
		node.sort_neighbors();

	message_stream << "Building priority queue..." << endl;
	queue = std::move(Queue(CmpEdge(), std::move(edge_vec)));
	message_stream << "#Edges: " << edges.size() << endl;

	return lambda;
}

void upgma() {
	const double max_dist = 10.0;
	message_stream << "Reading edges..." << endl;
	EdgeVec all_edges(config.query_file.c_str());
	message_stream << "Read " << all_edges.nodes() << " nodes, " << all_edges.size() << " edges." << endl;
	EdgeVec::iterator begin = all_edges.begin();

	EdgeList edges;
	vector<Node> nodes;
	for (int i = 0; i < (int)all_edges.nodes(); ++i)
		nodes.emplace_back(i, 1, i);
	Queue queue;
	double lambda = 0.0;
	int node_count = (int)nodes.size(), round = 0;
	do {
		lambda = load_edges(begin, all_edges.end(), edges, nodes, queue, lambda, max_dist);
		message_stream << "Clustering nodes..." << endl;
		message_stream << "#Edges: " << edges.size() << ", #Nodes: " << node_count << endl;
		while (!queue.empty()) {
			EdgePtr e = queue.top();
			queue.pop();
			while (!queue.empty() && !valid(queue.top(), nodes)) {
				EdgePtr f = queue.top();
				erase(f, edges);
				queue.pop();
			}
			if (!queue.empty() && !(*e <= *queue.top())) {
				//std::cerr << e->u << '\t' << queue.top()->l << '\t' << lambda << endl;
				queue.push(e);
				break;
			}
			if (valid(e, nodes) && e->u < max_dist) {
				merge_nodes(e->n1, e->n2, nodes, edges, queue, max_dist, lambda);
				--node_count;
				cout << nodes.back().parent << '\t' << e->n1 << '\t' << e->n2 << '\t' << std::setprecision(10) << e->u << '\t' << round << endl;
			}
			erase(e, edges);
			if (edges.size() % 10000 == 0)
				message_stream << "#Edges: " << edges.size() << ", #Nodes: " << node_count << endl;
		}
		message_stream << "#Edges: " << edges.size() << ", #Nodes: " << node_count << endl;
		++round;
	} while (lambda < max_dist);
}

}}}