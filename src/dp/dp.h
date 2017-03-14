/****
Copyright (c) 2016, University of Tuebingen, Benjamin Buchfink
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

#ifndef DP_H_
#define DP_H_

#include <utility>
#include <map>
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
	unsigned frame_, subject_, subject_pos_, query_pos_;
	Diagonal_segment ungapped;
	unsigned prefix_score;
};

template<typename _score>
void smith_waterman(const Letter *query, local_match &segment, _score gap_open, _score gap_extend, vector<char> &transcript_buf, const _score& = int());
int smith_waterman(const sequence &query, const sequence &subject, unsigned band, unsigned padding, int op, int ep);

int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned seed_len, unsigned &delta, unsigned &len);
int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned &delta, unsigned &len);
int xdrop_ungapped_right(const Letter *query, const Letter *subject, int &len);

struct Local {};
struct Global {};

void greedy_align(sequence query, const Long_score_profile &qp, sequence subject, vector<Seed_hit>::const_iterator begin, vector<Seed_hit>::const_iterator end, bool log, Hsp_data &out);
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
		for (int i = 0; i < buf.col_size_; ++i) {
			for (int j = 0; j < buf.data_.size() / buf.col_size_; ++j)
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
		link_idx(0)
	{}
	Diagonal_node(int query_pos, int subject_pos, int len, int score, int link_idx=0) :
		Diagonal_segment(query_pos, subject_pos, len, score),
		link_idx(link_idx)
	{}
	Diagonal_node(const Diagonal_segment &d) :
		Diagonal_segment(d),
		link_idx(0)
	{}
	unsigned link_idx;
};

struct Diag_graph
{
	struct Edge
	{
		Edge() :
			node_in(),
			exact(true)
		{
			prefix_score[0] = 0;
			prefix_score[1] = 0;
		}
		Edge(int prefix_score, int j, unsigned node_in, unsigned node_out, bool exact, unsigned state, int diff1, int diff2) :
			j(j),
			node_in(node_in),
			node_out(node_out),
			exact(exact),
			state(state),
			diff1(diff1),
			diff2(diff2)
		{
			this->prefix_score[0] = prefix_score;
			this->prefix_score[1] = 0;
		}
		/*operator int() const
		{
			return prefix_score;
		}
		bool operator<(const Edge &x) const
		{
			return prefix_score > x.prefix_score;
		}*/
		enum { n_pass = 2 };
		int prefix_score[n_pass], j, diff1, diff2;
		unsigned node_in, node_out, state;
		bool exact;
	};

	void init()
	{
		nodes.clear();
		edges.clear();
	}

	vector<Edge>::iterator add_edge(const Edge &edge)
	{
		for (vector<Diagonal_node>::iterator j = nodes.begin() + edge.node_in + 1; j < nodes.end(); ++j)
			++j->link_idx;
		return edges.insert(edges.begin() + nodes[edge.node_in].link_idx++, edge);
	}

	vector<Edge>::const_iterator get_edge(unsigned node, int j, int pass) const
	{
		const Diagonal_node &d = nodes[node];
		if (d.score == 0)
			return edges.begin() + d.link_idx - 1;
		if (edges.empty())
			return edges.end();
		int max_score = d.score;
		vector<Edge>::const_iterator max_edge = edges.end();
		for (vector<Edge>::const_iterator i = edges.begin() + d.link_idx - 1; i >= edges.begin() && i->node_in == node; --i)
			if (i->j < j && i->prefix_score[pass] > max_score) {
				max_edge = i;
				max_score = i->prefix_score[pass];
			}
		return max_edge;
	}

	int prefix_score(unsigned node, int j, int pass) const
	{
		const vector<Edge>::const_iterator i = get_edge(node, j, pass);
		return i == edges.end() ? nodes[node].score : std::max(nodes[node].score, i->prefix_score[pass]);
	}

	int prefix_score(unsigned node, int pass) const
	{
		return prefix_score(node, nodes[node].subject_end(), pass);
	}

	Diagonal_node& operator[](unsigned k)
	{
		return nodes[k];
	}

	const Diagonal_node& operator[](unsigned k) const
	{
		return nodes[k];
	}

	vector<Diagonal_node> nodes;
	vector<Edge> edges;
};

int needleman_wunsch(sequence query, sequence subject, int qbegin, int qend, int sbegin, int send, unsigned node, unsigned edge, Diag_graph &diags, bool log);

struct Band
{
	void init(int diags, int cols, bool zero)
	{
		diags_ = diags;
		cols_ = cols;
		data_.clear();
		data_.resize((size_t)diags*cols);
		if (zero)
			memset(data_.data(), 0, data_.size());
	}
	struct Iterator {
		Iterator(uint8_t *p, int diags) :
			p_(p),
			diags_(diags)
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
	uint8_t* data()
	{
		return data_.data();
	}
private:
	int diags_, cols_;
	vector<uint8_t> data_;
};

struct Diag_scores {
	enum {
		block_len = 16
	};
	void get_diag(int i, int j, int o, int j_begin, int j_end, vector<Diagonal_node> &diags, int cutoff, bool log);
	void get_diag2(int i, int j, int o, int j_begin, int j_end, vector<Diagonal_node> &diags, int cutoff, bool log);
	void scan_diags(int diag, sequence query, sequence subject, const Long_score_profile &qp, bool log, vector<Diagonal_node> &diags, std::multimap<int, unsigned> &window);
	void scan_vicinity(unsigned d_idx, unsigned e_idx, vector<Diagonal_node> &diags, bool log, uint8_t *sv_max);
	void set_zero(Band::Iterator &d, Band::Iterator d2, int begin, int end);
	void scan_ends(unsigned d_idx, vector<Diagonal_node> &diags, bool log, uint8_t *sv_max, int subject_len);
	void set_active(int o, int begin, int end);
	bool is_active(int o, int i) const;
	Band score_buf, local_max;
	vector<bool> active;
	int i_begin, j_begin, d_begin, d_end;
	static int min_diag_score, min_low_score;
};

void smith_waterman(sequence q, sequence s, Hsp_data &out);
void smith_waterman(sequence q, sequence s, const Diag_graph &diags);
int score_range(sequence query, sequence subject, int i, int j, int j_end);

#endif /* FLOATING_SW_H_ */
