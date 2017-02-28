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

template<typename _score>
void smith_waterman(const Letter *query, local_match &segment, _score gap_open, _score gap_extend, vector<char> &transcript_buf, const _score& = int());
int smith_waterman(const sequence &query, const sequence &subject, unsigned band, unsigned padding, int op, int ep);

int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned seed_len, unsigned &delta, unsigned &len);
int xdrop_ungapped(const Letter *query, const Letter *subject, unsigned &delta, unsigned &len);
int xdrop_ungapped_right(const Letter *query, const Letter *subject, int &len);

void greedy_align(sequence query, sequence subject, const vector<Diagonal_segment> &sh, bool log);
void greedy_align(sequence query, sequence subject, const Diagonal_segment &sh, bool log);
void greedy_align(sequence query, const Long_score_profile &qp, sequence subject, const Diagonal_segment &sh, bool log);
void greedy_align2(sequence query, const Long_score_profile &qp, sequence subject, const vector<Diagonal_segment> &sh, bool log, Hsp_data &out);
void greedy_align(sequence query, const Long_score_profile &qp, sequence subject, int qa, int sa, bool log);
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

struct Diagonal_node : public Diagonal_segment
{
	enum { n_path = 3, estimate, finished };
	struct Edge
	{
		Edge() :
			prefix_score(0),
			node(),
			exact(true)
		{}
		Edge(int prefix_score, int j, unsigned node, bool exact, unsigned state, int diff1, int diff2) :
			prefix_score(prefix_score),
			j(j),
			node(node),
			exact(exact),
			state(state),
			diff1(diff1),
			diff2(diff2)
		{}
		operator int() const
		{
			return prefix_score;
		}
		bool operator<(const Edge &x) const
		{
			return prefix_score > x.prefix_score;
		}
		int prefix_score, j, diff1, diff2;
		unsigned node, state;
		bool exact;
	};

	Diagonal_node() :
		Diagonal_segment()
	{}
	Diagonal_node(int query_pos, int subject_pos, int len, int score) :
		Diagonal_segment(query_pos, subject_pos, len, score)
	{}
	Diagonal_node(const Diagonal_segment &d) :
		Diagonal_segment(d)
	{}
	int prefix_score(int j) const
	{
		for (int k = 0; k < n_path; ++k)
			if (edges[k].j < j)
				return std::max(edges[k].prefix_score, score);
		return score;
	}
	int prefix_score() const
	{
		return std::max(edges[0].prefix_score, score);
	}
	int get_edge(int j) const
	{
		if (score == 0)
			return 0;
		int max_score = score, max_edge = -1;
		for (int k = 0; k < n_path; ++k)
			if (edges[k].j < j && edges[k].prefix_score > max_score) {
				max_edge = k;
				max_score = edges[k].prefix_score;
			}
		return max_edge;
	}
	Top_list<Edge, n_path> edges;
};

int needleman_wunsch(sequence query, sequence subject, int qbegin, int qend, int sbegin, int send, unsigned node, unsigned edge, vector<Diagonal_node> &diags, bool log);

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
		Iterator(const uint8_t *p, int diags) :
			p_(p),
			diags_(diags)
		{}
		uint8_t operator[](int i) const
		{
			return p_[i*diags_];
		}
	private:
		const int diags_;
		const uint8_t *p_;
	};
	Iterator diag(int o) const
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
	enum { block_len = 16 };
	void get_diag(int i, int j, int o, const sequence &query, const sequence &subject, vector<Diagonal_node> &diags);
	void scan_diags(const Diagonal_segment &diag, sequence query, sequence subject, const Long_score_profile &qp, bool log, vector<Diagonal_node> &diags);
	Band score_buf, local_max;
};

void smith_waterman(sequence q, sequence s, Hsp_data &out);
void smith_waterman(sequence q, sequence s, const vector<Diagonal_node> &diags);
int score_range(sequence query, sequence subject, int i, int j, int j_end);

#endif /* FLOATING_SW_H_ */
