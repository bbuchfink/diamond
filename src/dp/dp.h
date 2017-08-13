/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef DP_H_
#define DP_H_

#include <utility>
#include <map>
#include <list>
#include "../basic/match.h"
#include "../align/align.h"
#include "score_profile.h"

void init_cbs();

struct No_score_correction
{
	void operator()(int &score, int i, int query_anchor, int mult) const
	{}
};

struct Bias_correction : public vector<float>
{
	Bias_correction(const sequence &seq);
	void operator()(float &score, int i, int query_anchor, int mult) const
	{
		score += (*this)[query_anchor + i*mult];
	}
	int operator()(const Hsp_data &hsp) const;
	int operator()(const Diagonal_segment &d) const;
};

struct Seed_hit
{
	Seed_hit()
	{}
	Seed_hit(unsigned frame, unsigned subject, unsigned subject_pos, unsigned query_pos, const Diagonal_segment &ungapped) :
		frame_(frame),
		subject_(subject),
		subject_pos_(subject_pos),
		query_pos_(query_pos),
		ungapped(ungapped),
		prefix_score(ungapped.score)
	{ }
	int diagonal() const
	{
		return (int)query_pos_ - (int)subject_pos_;
	}
	bool operator<(const Seed_hit &rhs) const
	{
		return ungapped.score > rhs.ungapped.score;
	}
	static bool compare_pos(const Seed_hit &x, const Seed_hit &y)
	{
		return Diagonal_segment::cmp_subject_end(x.ungapped, y.ungapped);
	}
	static bool compare_diag(const Seed_hit &x, const Seed_hit &y)
	{
		return x.frame_ < y.frame_ || (x.frame_ == y.frame_ && (x.diagonal() < y.diagonal() || (x.diagonal() == y.diagonal() && x.ungapped.j < y.ungapped.j)));
	}
	struct Frame
	{
		unsigned operator()(const Seed_hit &x) const
		{
			return x.frame_;
		}
	};

	unsigned frame_, subject_, subject_pos_, query_pos_;
	Diagonal_segment ungapped;
	unsigned prefix_score;
};

struct Hsp_traits
{
	Hsp_traits(unsigned frame) :
		d_min(std::numeric_limits<int>::max()),
		d_max(std::numeric_limits<int>::min()),
		score(0),
		frame((int)frame)
	{}
	Hsp_traits(const interval &query_source_range):
		query_source_range(query_source_range)
	{}
	Hsp_traits(const Hsp_data &hsp)
	{}
	int partial_score(const Diagonal_segment &d) const
	{
		const double overlap = std::max(d.subject_range().overlap_factor(subject_range), d.query_range().overlap_factor(query_range));
		return int((1 - overlap)*d.score);
	}
	int partial_score(const Hsp_traits &x) const
	{
		const double overlap = std::max(x.subject_range.overlap_factor(subject_range), x.query_range.overlap_factor(query_range));
		return int((1 - overlap)*x.score);
	}
	bool disjoint(const Diagonal_segment &d) const
	{
		return intersect(query_range, d.query_range()).length() == 0 && intersect(subject_range, d.subject_range()).length() == 0;
	}
	bool disjoint(const Hsp_traits &x) const
	{
		return intersect(query_range, x.query_range).length() == 0 && intersect(subject_range, x.subject_range).length() == 0;
	}
	bool rel_disjoint(const Diagonal_segment &d) const
	{
		return intersect(query_range, d.query_range()).length() == 0 || intersect(subject_range, d.subject_range()).length() == 0;
	}
	bool rel_disjoint(const Hsp_traits &x) const
	{
		return intersect(query_range, x.query_range).length() == 0 || intersect(subject_range, x.subject_range).length() == 0;
	}
	bool collinear(const Hsp_traits &x) const
	{
		const int di = x.query_range.begin_ - query_range.begin_, dj = x.subject_range.begin_ - subject_range.begin_;
		return (di >= 0 && dj >= 0) || (di <= 0 && dj <= 0);
	}
	bool collinear(const Diagonal_segment &d) const
	{
		const int di = d.i - query_range.begin_, dj = d.j - subject_range.begin_;
		return (di >= 0 && dj >= 0) || (di <= 0 && dj <= 0);
	}
	static bool cmp_diag(const Hsp_traits &x, const Hsp_traits &y)
	{
		return x.frame < y.frame || (x.frame == y.frame && x.d_min < y.d_min);
	}
	struct Frame
	{
		unsigned operator()(const Hsp_traits &x) const
		{
			return x.frame;
		}
	};
	int d_min, d_max, score, frame;
	interval query_source_range, query_range, subject_range;
};

template<typename _score>
void smith_waterman(const Letter *query, local_match &segment, _score gap_open, _score gap_extend, vector<char> &transcript_buf, const _score& = int());
int smith_waterman(const sequence &query, const sequence &subject, unsigned band, unsigned padding, int op, int ep);

int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned seed_len, unsigned &delta, unsigned &len);
int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned &delta, unsigned &len);
int xdrop_ungapped_right(const Letter *query, const Letter *subject, int &len);
Diagonal_segment xdrop_ungapped(const sequence &query, const Bias_correction &query_bc, const sequence &subject, int qa, int sa);
Diagonal_segment xdrop_ungapped(const sequence &query, const sequence &subject, int qa, int sa);

struct Local {};
struct Global {};

int greedy_align(sequence query, const Long_score_profile &qp, const Bias_correction &query_bc, sequence subject, vector<Seed_hit>::const_iterator begin, vector<Seed_hit>::const_iterator end, bool log, std::list<Hsp_data> &hsps, std::list<Hsp_traits> &ts, unsigned frame);
int greedy_align(sequence query, const Long_score_profile &qp, const Bias_correction &query_bc, sequence subject, bool log, std::list<Hsp_data> &hsps, std::list<Hsp_traits>::const_iterator t_begin, std::list<Hsp_traits>::const_iterator t_end, std::list<Hsp_traits> &ts, int cutoff, unsigned frame);
int estimate_score(const Long_score_profile &qp, sequence s, int d, int d1, bool log);

template<typename _t>
struct Fixed_score_buffer
{

	inline void init(size_t col_size, size_t cols, _t init)
	{
		col_size_ = col_size;
		data_.clear();
		data_.reserve(col_size * cols);
		data_.resize(col_size);
		for (size_t i = 0; i<col_size; ++i)
			data_[i] = init;
	}
	
	std::pair<int, int> find(_t s) const
	{
		const int i = int(std::find(data_.begin(), data_.end(), s) - data_.begin());
		return std::pair<int, int>(int(i%col_size_), int(i / col_size_));
	}

	inline std::pair<_t*, _t*> get()
	{
		data_.resize(data_.size() + col_size_);
		_t* ptr = last();
		return std::pair<_t*, _t*>(ptr - col_size_, ptr);
	}

	inline _t* last()
	{
		return &*(data_.end() - col_size_);
	}

	const _t* column(int col) const
	{
		return &data_[col_size_*col];
	}

	_t operator()(int i, int j) const
	{
		return data_[j*col_size_ + i];
	}

	friend std::ostream& operator<<(std::ostream &s, const Fixed_score_buffer &buf)
	{
		s << '\t';
		for (int j = 0; j < int(buf.data_.size() / buf.col_size_); ++j)
			s << j << '\t';
		s << endl;
		for (int i = 0; i < int(buf.col_size_); ++i) {
			s << i << '\t';
			for (int j = 0; j < int(buf.data_.size() / buf.col_size_); ++j)
				s << buf(i, j) << '\t';
			s << endl;
		}
		return s;
	}

private:
	vector<_t> data_;
	size_t col_size_;

};

template<typename _score, typename _mode>
const Fixed_score_buffer<_score>& needleman_wunsch(sequence query, sequence subject, int &max_score, const _mode&, const _score&);

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
	Diagonal_node(int query_pos, int subject_pos, int len, int score, int link_idx=-1) :
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

	enum { end = 0xffffffffffffffffllu };

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

	void load(vector<Seed_hit>::const_iterator begin, vector<Seed_hit>::const_iterator end);
	void sort();
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
			if (i->j <= j && i->prefix_score > max_score) {
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

int needleman_wunsch(sequence query, sequence subject, int qbegin, int qend, int sbegin, int send, unsigned node, unsigned edge, Diag_graph &diags, bool log);

struct Band
{
	void init(int diags, int cols)
	{
		diags_ = diags;
		cols_ = cols;
		data_.clear();
		data_.resize((size_t)diags*cols);
	}
	struct Iterator {
		Iterator(uint8_t *p, int diags) :
			diags_(diags),
			p_(p)			
		{}
		uint8_t& operator[](int i)
		{
			return p_[i*diags_];
		}
	private:
		const int diags_;
		uint8_t *p_;
	};
	Iterator diag(int o)
	{
		return Iterator(&data_[o], diags_);
	}
	int cols() const
	{
		return cols_;
	}
	int diags() const
	{
		return diags_;
	}
	uint8_t* data()
	{
		return data_.data();
	}
	bool check(uint8_t *ptr) const
	{
		return ptr >= data_.data() && ptr <= data_.data() + data_.size();
	}
private:
	int diags_, cols_;
	vector<uint8_t> data_;
};

struct Diag_scores {
	enum {
		block_len = 16
	};
	int dj0(int d) const
	{
		return std::max(-d, 0);
	}
	int dj1(int d) const
	{
		return std::min(qlen - d, slen);
	}
	void get_diag(int i, int j, int o, int j_begin, int j_end, vector<Diagonal_node> &diags, int cutoff, bool log, const Bias_correction &query_bc);
	void scan_diags(int d_begin, int d_end, sequence query, sequence subject, const Long_score_profile &qp, const Bias_correction &query_bc, bool log, vector<Diagonal_node> &diags, bool fast);
	Band score_buf, local_max;
	vector<uint8_t> sv_max;
	vector<bool> active;
	int i_begin, j_begin, d_begin, d_end, qlen, slen;
	bool fast;
	static int min_diag_score, min_low_score;
};

void smith_waterman(sequence q, sequence s, Hsp_data &out);
void smith_waterman(sequence q, sequence s, const Diag_graph &diags);
int score_range(sequence query, sequence subject, int i, int j, int j_end);

void swipe(const sequence &query, vector<sequence>::const_iterator subject_begin, vector<sequence>::const_iterator subject_end, vector<int>::iterator out);
void banded_sw(const sequence &query, const sequence &subject, int d_begin, int d_end, int j_begin, int j_end, Hsp_data &out);

#endif /* FLOATING_SW_H_ */
