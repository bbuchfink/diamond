#pragma once

#include <vector>
#include <algorithm>

/* Suppose there are N=2^(K+1)-1 sorted numbers in an array a[]. They
 * implicitly form a complete binary tree of height K+1. We consider leaves to
 * be at level 0. The binary tree has the following properties:
 *
 * 1. The lowest k-1 bits of nodes at level k are all 1. The k-th bit is 0.
 *    The first node at level k is indexed by 2^k-1. The root of the tree is
 *    indexed by 2^K-1.
 *
 * 2. For a node x at level k, its left child is x-2^(k-1) and the right child
 *    is x+2^(k-1).
 *
 * 3. For a node x at level k, it is a left child if its (k+1)-th bit is 0. Its
 *    parent node is x+2^k. Similarly, if the (k+1)-th bit is 1, x is a right
 *    child and its parent is x-2^k.
 *
 * 4. For a node x at level k, there are 2^(k+1)-1 nodes in the subtree
 *    descending from x, including x. The left-most leaf is x&~(2^k-1) (masking
 *    the lowest k bits to 0).
 *
 * When numbers can't fill a complete binary tree, the parent of a node may not
 * be present in the array. The implementation here still mimics a complete
 * tree, though getting the special casing right is a little complex. There may
 * be alternative solutions.
 *
 * As a sorted array can be considered as a binary search tree, we can
 * implement an interval tree on top of the idea. We only need to record, for
 * each node, the maximum value in the subtree descending from the node.
 */
template<typename S, typename T> // "S" is a scalar type; "T" is the type of data associated with each interval
class IITree {
	struct StackCell {
		size_t x; // node
		int k, w; // k: level; w: 0 if left child hasn't been processed
		StackCell() {};
		StackCell(int k_, size_t x_, int w_) : x(x_), k(k_), w(w_) {};
	};
	struct Interval {
		S st, en, max;
		T data;
		Interval(const S &s, const S &e, const T &d) : st(s), en(e), max(e), data(d) {};
	};
	struct IntervalLess {
		bool operator()(const Interval &a, const Interval &b) const { return a.st < b.st; }
	};
	std::vector<Interval> a;
	int max_level;
	int index_core(std::vector<Interval> &a) {
		size_t i, last_i; // last_i points to the rightmost node in the tree
		S last; // last is the max value at node last_i
		int k;
		if (a.size() == 0) return -1;
		for (i = 0; i < a.size(); i += 2) last_i = i, last = a[i].max = a[i].en; // leaves (i.e. at level 0)
		for (k = 1; 1LL<<k <= a.size(); ++k) { // process internal nodes in the bottom-up order
			size_t x = 1LL<<(k-1), i0 = (x<<1) - 1, step = x<<2; // i0 is the first node
			for (i = i0; i < a.size(); i += step) { // traverse all nodes at level k
				S el = a[i - x].max;                          // max value of the left child
				S er = i + x < a.size()? a[i + x].max : last; // of the right child
				S e = a[i].en;
				e = e > el? e : el;
				e = e > er? e : er;
				a[i].max = e; // set the max value for node i
			}
			last_i = last_i>>k&1? last_i - x : last_i + x; // last_i now points to the parent of the original last_i
			if (last_i < a.size() && a[last_i].max > last) // update last accordingly
				last = a[last_i].max;
		}
		return k - 1;
	}
public:
	void add(const S &s, const S &e, const T &d) { a.push_back(Interval(s, e, d)); }
	void index(void) {
		std::sort(a.begin(), a.end(), IntervalLess());
		max_level = index_core(a);
	}
	bool overlap(const S &st, const S &en, std::vector<size_t> &out) const {
		int t = 0;
		StackCell stack[64];
		out.clear();
		if (max_level < 0) return false;
		stack[t++] = StackCell(max_level, (1LL<<max_level) - 1, 0); // push the root; this is a top down traversal
		while (t) { // the following guarantees that numbers in out[] are always sorted
			StackCell z = stack[--t];
			if (z.k <= 3) { // we are in a small subtree; traverse every node in this subtree
				size_t i, i0 = z.x >> z.k << z.k, i1 = i0 + (1LL<<(z.k+1)) - 1;
				if (i1 >= a.size()) i1 = a.size();
				for (i = i0; i < i1 && a[i].st < en; ++i)
					if (st < a[i].en) // if overlap, append to out[]
						out.push_back(i);
			} else if (z.w == 0) { // if left child not processed
				size_t y = z.x - (1LL<<(z.k-1)); // the left child of z.x; NB: y may be out of range (i.e. y>=a.size())
				stack[t++] = StackCell(z.k, z.x, 1); // re-add node z.x, but mark the left child having been processed
				if (y >= a.size() || a[y].max > st) // push the left child if y is out of range or may overlap with the query
					stack[t++] = StackCell(z.k - 1, y, 0);
			} else if (z.x < a.size() && a[z.x].st < en) { // need to push the right child
				if (st < a[z.x].en) out.push_back(z.x); // test if z.x overlaps the query; if yes, append to out[]
				stack[t++] = StackCell(z.k - 1, z.x + (1LL<<(z.k-1)), 0); // push the right child
			}
		}
		return out.size() > 0? true : false;
	}
	size_t size(void) const { return a.size(); }
	const S &start(size_t i) const { return a[i].st; }
	const S &end(size_t i) const { return a[i].en; }
	const T &data(size_t i) const { return a[i].data; }
};
