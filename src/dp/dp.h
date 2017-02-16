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

	struct Edge
	{
		Edge() :
			prefix_score(0),
			node(),
			exact(true)
		{}
		Edge(int prefix_score, int j, unsigned node, bool exact) :
			prefix_score(prefix_score),
			j(j),
			node(node),
			exact(exact)
		{}
		operator int() const
		{
			return prefix_score;
		}
		int prefix_score, j;
		unsigned node;
		bool exact;
	};

	Diagonal_node() :
		Diagonal_segment(),
		diff(std::numeric_limits<int>::min())
	{}
	Diagonal_node(int query_pos, int subject_pos, int len, int score) :
		Diagonal_segment(query_pos, subject_pos, len, score),
		diff(std::numeric_limits<int>::min())
	{}
	Diagonal_node(const Diagonal_segment &d) :
		Diagonal_segment(d),
		diff(std::numeric_limits<int>::min())
	{}
	enum { n_path = 2 };
	Top_list<Edge, n_path> edges;
	int diff;
};

int needleman_wunsch(sequence query, sequence subject, int qbegin, int qend, int sbegin, int send, unsigned node, unsigned edge, vector<Diagonal_node> &diags, bool log);

#endif /* FLOATING_SW_H_ */
