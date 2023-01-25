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
#include "double_buffer.h"
#include "../../util/util.h"
#include "traceback.h"
#include "../output/output_format.h"

using std::vector;
using std::pair;

struct Local {};
struct Global {};

template<typename _score, typename _mode>
_score saturate(_score x)
{
	return x;
}

template<>
int saturate<int, Local>(int x)
{
	return std::max(x, 0);
}

template<typename _score, typename _mode>
void set_max_score(_score s, _score& max_score)
{
}

template<>
void set_max_score<int, Local>(int s, int& max_score)
{
	max_score = std::max(max_score, s);
}

template<typename _t>
struct Fixed_score_buffer
{

	inline void init(size_t col_size, size_t cols, _t init)
	{
		col_size_ = col_size;
		data_.clear();
		data_.reserve(col_size * cols);
		data_.resize(col_size);
		for (size_t i = 0; i < col_size; ++i)
			data_[i] = init;
	}

	std::pair<int, int> find(_t s) const
	{
		const int i = int(std::find(data_.begin(), data_.end(), s) - data_.begin());
		return std::pair<int, int>(int(i % col_size_), int(i / col_size_));
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
		return &data_[col_size_ * col];
	}

	_t operator()(int i, int j) const
	{
		return data_[j * col_size_ + i];
	}

	friend std::ostream& operator<<(std::ostream& s, const Fixed_score_buffer& buf)
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
	std::vector<_t> data_;
	size_t col_size_;

};

template<typename _score, typename _mode>
struct Dp_matrix
{

	struct Column_iterator
	{

		inline Column_iterator(const pair<_score*, _score*>& score, _score* hgap, int query_len, int col) :
			score_(score),
			hgap_(hgap),
			end_(score_.second + query_len + 1),
			i_(0)
		{
			*score_.first = saturate<_score, _mode>(col == 0 ? 0 : -score_matrix.gap_open() - col * score_matrix.gap_extend());
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

		inline _score& score()
		{
			return *score_.second;
		}

		inline _score diag() const
		{
			return *score_.first;
		}

		inline _score& hgap()
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
		pair<_score*, _score*> score_;
		_score* hgap_;
		const _score* const end_;
		int i_;

	};

	inline Column_iterator column(int j)
	{
		return Column_iterator(score_.get(), hgap_.data(), query_len_, j);
	}

	inline Dp_matrix(int query_len, int subject_len) :
		query_len_(query_len)
	{
		score_.init(query_len + 1, subject_len + 1, 0);
		hgap_.clear();
		hgap_.insert(hgap_.end(), query_len, std::numeric_limits<int>::min() + score_matrix.gap_extend());
		int* score = score_.last();
		int g = -score_matrix.gap_open() - score_matrix.gap_extend();
		for (int i = 1; i <= query_len; ++i)
			score[i] = saturate<_score, _mode>(g--);
	}

	const Fixed_score_buffer<_score>& score_buffer() const
	{
		return score_;
	}

private:

	const int query_len_;
	static thread_local Fixed_score_buffer<_score> score_;
	static thread_local vector<_score> hgap_;

};

template<typename _score, typename _mode> thread_local Fixed_score_buffer<_score> Dp_matrix<_score, _mode>::score_;
template<typename _score, typename _mode> thread_local vector<_score> Dp_matrix<_score, _mode>::hgap_;

template<typename _score, typename _mode>
const Fixed_score_buffer<_score>& needleman_wunsch(Sequence query, Sequence subject, int& max_score, const _mode&, const _score&)
{
	using std::max;
	const int gap_open = score_matrix.gap_open() + score_matrix.gap_extend(), gap_extend = score_matrix.gap_extend();
	int m = 0;

	Dp_matrix<_score, _mode> mtx((unsigned)query.length(), (unsigned)subject.length());

	for (int j = 0; j < (int)subject.length(); ++j) {
		typename Dp_matrix<_score, _mode>::Column_iterator it = mtx.column(j);
		_score vgap = std::numeric_limits<int>::min() + gap_extend;
		for (; it.valid(); ++it) {
			const _score match_score = score_matrix(subject[j], query[it.row()]);
			const _score s = saturate<_score, _mode>(max(max(it.diag() + match_score, vgap), it.hgap()));
			const _score open = s - gap_open;
			vgap = max(vgap - gap_extend, open);
			it.hgap() = max(it.hgap() - gap_extend, open);
			it.score() = s;
			set_max_score<_score, _mode>(s, m);
		}
	}

	max_score = m;
	return mtx.score_buffer();
}

template const Fixed_score_buffer<int>& needleman_wunsch<int, Local>(Sequence query, Sequence subject, int& max_score, const Local&, const int&);

void smith_waterman(Sequence q, Sequence s, Hsp& out)
{
	int max_score;
	const Fixed_score_buffer<int>& dp = needleman_wunsch(q, s, max_score, Local(), int());
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