/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#include <vector>
#include "../dp.h"
#include "traceback.h"

using std::vector;
using std::pair;

struct Local {};
struct Global {};

template<typename Score, typename Mode>
Score saturate(Score x)
{
	return x;
}

template<>
int saturate<int, Local>(int x)
{
	return std::max(x, 0);
}

template<typename Score, typename Mode>
void set_max_score(Score s, Score& max_score)
{
}

template<>
void set_max_score<int, Local>(int s, int& max_score)
{
	max_score = std::max(max_score, s);
}

template<typename T>
struct FixedScoreBuffer
{

	inline void init(size_t col_size, size_t cols, T init)
	{
		col_size_ = col_size;
		data_.clear();
		data_.reserve(col_size * cols);
		data_.resize(col_size);
		for (size_t i = 0; i < col_size; ++i)
			data_[i] = init;
	}

	std::pair<int, int> find(T s) const
	{
		const int i = int(std::find(data_.begin(), data_.end(), s) - data_.begin());
		return std::pair<int, int>(int(i % col_size_), int(i / col_size_));
	}

	inline std::pair<T*, T*> get()
	{
		data_.resize(data_.size() + col_size_);
		T* ptr = last();
		return std::pair<T*, T*>(ptr - col_size_, ptr);
	}

	inline T* last()
	{
		return &*(data_.end() - col_size_);
	}

	const T* column(int col) const
	{
		return &data_[col_size_ * col];
	}

	T operator()(int i, int j) const
	{
		return data_[j * col_size_ + i];
	}

	friend std::ostream& operator<<(std::ostream& s, const FixedScoreBuffer& buf)
	{
		s << '\t';
		for (int j = 0; j < int(buf.data_.size() / buf.col_size_); ++j)
			s << j << '\t';
		s << std::endl;
		for (int i = 0; i < int(buf.col_size_); ++i) {
			s << i << '\t';
			for (int j = 0; j < int(buf.data_.size() / buf.col_size_); ++j)
				s << buf(i, j) << '\t';
			s << std::endl;
		}
		return s;
	}

private:
	std::vector<T> data_;
	size_t col_size_;

};

template<typename Score, typename Mode>
struct DpMatrix
{

	struct Column_iterator
	{

		inline Column_iterator(const pair<Score*, Score*>& score, Score* hgap, int query_len, int col) :
			score_(score),
			hgap_(hgap),
			end_(score_.second + query_len + 1),
			i_(0)
		{
			*score_.first = saturate<Score, Mode>(col == 0 ? 0 : -score_matrix.gap_open() - col * score_matrix.gap_extend());
			++score_.second;
		}

		inline int row() const
		{
			return i_;
		}

		inline bool valid() const
		{
			return score_.second < end_;
		}

		inline Score& score()
		{
			return *score_.second;
		}

		inline Score diag() const
		{
			return *score_.first;
		}

		inline Score& hgap()
		{
			return *hgap_;
		}

		inline void operator++()
		{
			++i_;
			++score_.first;
			++score_.second;
			++hgap_;
		}

	private:
		pair<Score*, Score*> score_;
		Score* hgap_;
		const Score* const end_;
		int i_;

	};

	inline Column_iterator column(int j)
	{
		return Column_iterator(score_.get(), hgap_.data(), query_len_, j);
	}

	inline DpMatrix(int query_len, int subject_len) :
		query_len_(query_len)
	{
		score_.init(query_len + 1, subject_len + 1, 0);
		hgap_.clear();
		hgap_.insert(hgap_.end(), query_len, std::numeric_limits<int>::min() + score_matrix.gap_extend());
		int* score = score_.last();
		int g = -score_matrix.gap_open() - score_matrix.gap_extend();
		for (int i = 1; i <= query_len; ++i)
			score[i] = saturate<Score, Mode>(g--);
	}

	const FixedScoreBuffer<Score>& score_buffer() const
	{
		return score_;
	}

private:

	const int query_len_;
	FixedScoreBuffer<Score> score_;
	vector<Score> hgap_;

};

template<typename Score, typename Mode>
FixedScoreBuffer<Score> needleman_wunsch(Sequence query, Sequence subject, int& max_score, const Mode&, const Score&)
{
	using std::max;
	const int gap_open = score_matrix.gap_open() + score_matrix.gap_extend(), gap_extend = score_matrix.gap_extend();
	int m = 0;

	DpMatrix<Score, Mode> mtx((unsigned)query.length(), (unsigned)subject.length());

	for (int j = 0; j < (int)subject.length(); ++j) {
		typename DpMatrix<Score, Mode>::Column_iterator it = mtx.column(j);
		Score vgap = std::numeric_limits<int>::min() + gap_extend;
		for (; it.valid(); ++it) {
			const Score match_score = score_matrix(subject[j], query[it.row()]);
			const Score s = saturate<Score, Mode>(max(max(it.diag() + match_score, vgap), it.hgap()));
			const Score open = s - gap_open;
			vgap = max(vgap - gap_extend, open);
			it.hgap() = max(it.hgap() - gap_extend, open);
			it.score() = s;
			set_max_score<Score, Mode>(s, m);
		}
	}

	max_score = m;
	return mtx.score_buffer();
}

template FixedScoreBuffer<int> needleman_wunsch<int, Local>(Sequence query, Sequence subject, int& max_score, const Local&, const int&);

void smith_waterman(Sequence q, Sequence s, Hsp& out)
{
	int max_score;
	const FixedScoreBuffer<int> dp = needleman_wunsch(q, s, max_score, Local(), int());
	pair<int, int> max_pos = dp.find(max_score);

	const int gap_open = score_matrix.gap_open(), gap_extend = score_matrix.gap_extend();
	int l, i = max_pos.first, j = max_pos.second, score;
	out.clear();
	out.score = dp(i, j);
	out.query_range.end_ = i;
	out.subject_range.end_ = j;

	while ((score = dp(i, j)) > 0) {
		const int match_score = score_matrix(q[i - 1], s[j - 1]);
		if (score == match_score + dp(i - 1, j - 1)) {
			if (q[i - 1] == s[j - 1]) {
				out.transcript.push_back(op_match);
				++out.identities;
			}
			else {
				out.transcript.push_back(op_substitution, s[j - 1]);
			}
			--i;
			--j;
			++out.length;
		}
		else if (have_hgap(dp, i, j, gap_open, gap_extend, l)) {
			for (; l > 0; l--) {
				out.transcript.push_back(op_deletion, s[--j]);
				++out.length;
			}
		}
		else if (have_vgap(dp, i, j, gap_open, gap_extend, l)) {
			out.transcript.push_back(op_insertion, (unsigned)l);
			out.length += l;
			i -= l;
		}
		else
			throw std::runtime_error("Traceback error.");
	}

	out.query_range.begin_ = i;
	out.subject_range.begin_ = j;
	out.query_source_range = out.query_range;
	out.transcript.reverse();
	out.transcript.push_terminator();
}