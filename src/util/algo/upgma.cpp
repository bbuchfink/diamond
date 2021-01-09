#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <limits.h>
#include "../../lib/MemoryPool/MemoryPool.h"
#include "../io/text_input_file.h"
#include "../../basic/config.h"
#include "../string/tokenizer.h"
#include "../log_stream.h"

using std::string;
using std::list;
using std::endl;
using std::map;
using std::vector;
using std::ostream;

namespace Util { namespace Algo { namespace UPGMA {

struct Edge {
	Edge(int n1, int n2, double d):
		n1(n1),
		n2(n2),
		deleted(0),
		d(d)
	{}
	bool operator<(const Edge &e) const {
		return d < e.d;
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
		s << e.n1 << '\t' << e.n2 << '\t' << e.d;
		return s;
	}
	int n1, n2, deleted;
	double d;
};

typedef list<Edge, MemoryPool<Edge>> EdgeList;
typedef EdgeList::iterator EdgePtr;

struct CmpEdge {
	bool operator()(const EdgePtr& e, const EdgePtr &f) const {
		return e->d > f->d;
	}
};

typedef std::priority_queue<EdgePtr, vector<EdgePtr>, CmpEdge> Queue;

struct Node {
	Node(int size, int parent):
		size(size),
		parent(parent)
	{}
	void sort_neighbors() {
		std::sort(neighbors.begin(), neighbors.end(), CmpNeighbor{ parent });
	}
	void set_parent(int parent, EdgeList &edges) {
		this->parent = parent;
		for (EdgePtr e : neighbors) {
			++e->deleted;
			if (e->deleted == 3)
				edges.erase(e);
		}
		neighbors.clear();
		neighbors.shrink_to_fit();
	}
	struct CmpNeighbor {
		int me;
		bool operator()(const EdgePtr &e, const EdgePtr &f) const { 
			return e->target(me) < f->target(me);
		}
	};
	int size;
	int parent;
	vector<EdgePtr> neighbors;
};

void merge_nodes(int n1,
	int n2,
	vector<Node> &nodes,
	EdgeList &edges,
	Queue &queue,
	double max_dist) {	
	const int union_idx = (int)nodes.size();
	nodes.emplace_back(nodes[n1].size + nodes[n2].size, union_idx);
	Node &node1 = nodes[n1], &node2 = nodes[n2], &union_node = nodes.back();
	
	vector<EdgePtr>::iterator i = node1.neighbors.begin(), j = node2.neighbors.begin();
	while (i < node1.neighbors.end() || j < node2.neighbors.end()) {
		int it = i < node1.neighbors.end() ? (*i)->target(n1) : INT_MAX, jt = j < node2.neighbors.end() ? (*j)->target(n2) : INT_MAX;
		double d1, d2;
		if (it == jt) {
			d1 = (*i)->d;
			d2 = (*j)->d;
			++i;
			++j;
		}
		else if (it < jt) {
			d1 = (*i)->d;
			d2 = max_dist;
			++i;
		}
		else {
			d1 = max_dist;
			d2 = (*j)->d;
			it = jt;
			++j;
		}
		if (nodes[it].parent != it || it == n1 || it == n2)
			continue;
		const EdgePtr e = edges.emplace(edges.end(), it, union_idx, (node1.size * d1 + node2.size * d2) / (node1.size + node2.size));
		queue.push(e);
		union_node.neighbors.push_back(e);
		nodes[it].neighbors.push_back(e);
	}
	node1.set_parent(union_idx, edges);
	node2.set_parent(union_idx, edges);
}

void upgma(EdgeList &edges, size_t node_count) {
	const double max_dist = 10.0;
	vector<EdgePtr> edge_vec;
	vector<Node> nodes;

	message_stream << "Building node vector..." << endl;
	nodes.reserve(node_count);
	for (int i = 0; i<int(node_count); ++i)
		nodes.emplace_back(1, i);

	message_stream << "Building edge vector and neighborhood..." << endl;
	edge_vec.reserve(edges.size());
	for (EdgePtr i = edges.begin(); i != edges.end(); ++i) {
		edge_vec.push_back(i);
		nodes[i->n1].neighbors.push_back(i);
		nodes[i->n2].neighbors.push_back(i);
	}

	message_stream << "Sorting neighborhoods..." << endl;
	for (Node &node : nodes)
		node.sort_neighbors();
	
	message_stream << "Building priority queue..." << endl;
	Queue queue(CmpEdge(), std::move(edge_vec));

	message_stream << "Clustering nodes..." << endl;
	while (!queue.empty()) {
		EdgePtr e = queue.top();
		queue.pop();
		if (nodes[e->n1].parent == e->n1 && nodes[e->n2].parent == e->n2 && e->d < max_dist) {
			merge_nodes(e->n1, e->n2, nodes, edges, queue, max_dist);
			--node_count;
		}
		++e->deleted;
		if (e->deleted == 3)
			edges.erase(e);
		if (edges.size() % 10000 == 0)
			message_stream << "#Edges: " << edges.size() << ", #Nodes: " << node_count << endl;
	}
}

void upgma() {
	TextInputFile in(config.query_file.front());
	string query, target;
	double evalue;
	EdgeList edges;	
	map<string, int> acc2idx;

	message_stream << "Reading edges..." << endl;
	while (in.getline(), !in.eof()) {
		String::Tokenizer t(in.line, "\t");
		t >> query >> target >> evalue;
		if (acc2idx.find(query) == acc2idx.end())
			acc2idx[query] = acc2idx.size();
		if (acc2idx.find(target) == acc2idx.end())
			acc2idx[target] = acc2idx.size();
		const int query_idx = acc2idx[query], target_idx = acc2idx[target];
		if (query_idx < target_idx)
			edges.emplace_back(query_idx, target_idx, evalue);
		if (edges.size() % 10000 == 0 && !edges.empty())
			message_stream << "#Edges read: " << edges.size() << endl;
	}
	message_stream << "#Edges: " << edges.size() << ", #Nodes: " << acc2idx.size() << endl;
	in.close();
	upgma(edges, acc2idx.size());
}

}}}