#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <limits.h>
#include "../io/text_input_file.h"
#include "../../basic/config.h"
#include "../string/tokenizer.h"
#include "../log_stream.h"
#include "../../lib/MemoryPool/MemoryPool.h"

using std::string;
using std::list;
using std::endl;
using std::map;
using std::vector;
using std::ostream;
using std::pair;

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
		return e->l > f->l;
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
	EdgePtr find_edge(int target, EdgePtr end) const {
		for (EdgePtr e : neighbors)
			if (e->n1 == target || e->n2 == target)
				return e;
		return end;
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
	double max_dist,
	double lambda) {	
	const int union_idx = (int)nodes.size();
	nodes.emplace_back(nodes[n1].size + nodes[n2].size, union_idx);
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

int node_idx(const string &acc, vector<Node> &nodes, map<string, int> &acc2idx, int &node_count) {
	if (acc2idx.find(acc) == acc2idx.end()) {
		const int i = (int)acc2idx.size();
		acc2idx[acc] = i;
		nodes.emplace_back(1, i);
		++node_count;
		return i;
	}
	else
		return acc2idx[acc];
}

double load_edges(TextInputFile &in, EdgeList &edges, vector<Node> &nodes, map<string, int> &acc2idx, Queue &queue, double lambda, double max_dist, int &node_count) {
	string query, target;
	double evalue = max_dist;
	message_stream << "Reading edges..." << endl;
	while (in.getline(), (!in.eof() && edges.size() < config.upgma_edge_limit)) {
		String::Tokenizer t(in.line, "\t");
		t >> query >> target >> evalue;
		const int query_idx = node_idx(query, nodes, acc2idx, node_count), target_idx = node_idx(target, nodes, acc2idx, node_count);

		if (query_idx < target_idx) {
			const int i = parent(query_idx, nodes), j = parent(target_idx, nodes);
			if (i == query_idx && j == target_idx)
				edges.emplace_back(i, j, 1, evalue);
			else {
				const EdgePtr e = nodes[i].find_edge(j, edges.end());
				const double max_edges = (double)nodes[i].size*(double)nodes[j].size;
				if (e == edges.end())
					edges.emplace(edges.end(), i, j, 1, evalue); // ->set_bounds(lambda, max_dist, max_edges);
				else {
					++e->count;
					e->s += evalue;
					//e->set_bounds(lambda, max_dist, max_edges);
				}
			}
		}

		if (edges.size() % 10000 == 0 && !edges.empty())
			message_stream << "#Edges: " << edges.size() << endl;
	}

	lambda = in.eof() ? max_dist : evalue;
	
	message_stream << "Clearing old edges..." << endl;
	for (EdgeList::iterator i = edges.begin(); i != edges.end();)
		if (nodes[i->n1].parent != i->n1 || nodes[i->n2].parent != i->n2)
			i = edges.erase(i);
		else {
			i->deleted = 0;
			i->set_bounds(lambda, max_dist, (double)nodes[i->n1].size * (double)nodes[i->n2].size);
			++i;
		}

	message_stream << "Clearing neighborhoods..." << endl;
	for (Node &node : nodes)
		node.neighbors.clear();

	message_stream << "Building edge vector and neighborhood..." << endl;
	vector<EdgePtr> edge_vec;
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
	queue = std::move(Queue(CmpEdge(), std::move(edge_vec)));
	message_stream << "#Edges: " << edges.size() << endl;

	return lambda;
}

void upgma() {
	const double max_dist = 10.0;
	TextInputFile in(config.query_file);
	EdgeList edges;
	vector<Node> nodes;
	map<string, int> acc2idx;
	Queue queue;
	double lambda = 0.0;
	int node_count = 0;
	do {
		lambda = load_edges(in, edges, nodes, acc2idx, queue, lambda, max_dist, node_count);
		message_stream << "Clustering nodes..." << endl;
		message_stream << "#Edges: " << edges.size() << ", #Nodes: " << node_count << endl;
		while (!queue.empty()) {
			EdgePtr e = queue.top();
			queue.pop();
			if (!queue.empty() && !(*e <= *queue.top())) {
				std::cout << e->u << '\t' << queue.top()->l << '\t' << lambda << endl;
				queue.push(e);
				break;
			}
			if (nodes[e->n1].parent == e->n1 && nodes[e->n2].parent == e->n2 && e->u < max_dist) {
				merge_nodes(e->n1, e->n2, nodes, edges, queue, max_dist, lambda);
				--node_count;
			}
			++e->deleted;
			if (e->deleted == 3)
				edges.erase(e);
			if (edges.size() % 10000 == 0)
				message_stream << "#Edges: " << edges.size() << ", #Nodes: " << node_count << endl;
		}
		message_stream << "#Edges: " << edges.size() << ", #Nodes: " << node_count << endl;
	} while (lambda < max_dist);
	in.close();	
}

}}}