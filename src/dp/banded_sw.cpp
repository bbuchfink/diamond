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

#include "dp.h"

using std::pair;

struct Banded_traceback_matrix
{
	Banded_traceback_matrix(const vector<int> &data, int band, int i0) :
		data_(data),
		band_(band),
		i0_(i0)
	{ }
	int operator()(int i, int j) const
	{
		return data_[(j + 1)*band_ + (i - (i0_ + j))];
	}
	bool in_band(int i, int j) const
	{
		return i >= 0 && j >= 0 && i >= i0_ + j && i < i0_ + j + band_;
	}
	void print(int qlen, int slen) const
	{
		printf("\n    ");
		for (int j = 0; j <= slen; ++j)
			printf("%4i", j - 1);
		printf("\n");
		for (int i = 0; i <= qlen; ++i) {
			printf("%4i", i - 1);
			for (int j = 0; j <= slen; ++j)
				printf("%4i", in_band(i - 1, j - 1) ? this->operator()(i - 1, j - 1) : 0);
			printf("\n");
		}
	}
	struct Column_iterator
	{
		Column_iterator(const int *ptr, const int *end):
			ptr_(ptr),
			end_(end)
		{}
		bool good()
		{
			return ptr_ >= end_;
		}
		void operator--()
		{
			--ptr_;
		}
		int operator*() const
		{
			return *ptr_;
		}
		const int *ptr_, *end_;
	};
	Column_iterator column(int i, int j) const
	{
		const int i0 = i0_ + j;
		return Column_iterator(&data_[(j + 1)*band_ + (i - i0)], &data_[(j + 1)*band_ + std::max(i0, 0) - i0]);
	}
	struct Row_iterator
	{
		Row_iterator(const int *ptr, const int *end, int band) :
			ptr_(ptr),
			end_(end),
			band_(band-1)
		{}
		bool good()
		{
			return ptr_ >= end_;
		}
		void operator--()
		{
			ptr_ -= band_;
		}
		int operator*() const
		{
			//cout << "*h=" << *ptr_ << endl;
			//printf("ptr=%llx end=%llx\n", ptr_, end_);
			return *ptr_;
		}
		const int *ptr_, *end_, band_;
	};
	Row_iterator row(int i, int j) const
	{
		const int i0 = i0_ + j;
		//cout << "j_end=" << i - i0_ - band_ << endl;
		const int* p = &data_[(j + 1)*band_ + (i - i0)];
		return Row_iterator(p, std::max(p - (j - (i - i0_ - band_) - 1)*(band_ - 1), &data_[band_]), band_);
	}
private:
	const vector<int> &data_;
	const int band_, i0_;
};

int have_gap(const Banded_traceback_matrix &dp,
	int i,
	int j,
	int &l)
{
	const int score = dp(i, j);
	l = 1;
	const int ge = score_matrix.gap_extend();
	int g = score_matrix.gap_open() + ge;
	Banded_traceback_matrix::Column_iterator v(dp.column(i - 1, j));
	Banded_traceback_matrix::Row_iterator h(dp.row(i, j - 1));
	while (v.good() && h.good()) {
		if (score == *v - g)
			return 0;
		else if (score == *h - g)
			return 1;
		--h;
		--v;
		++l;
		g += ge;
	}
	while (v.good()) {
		if (score == *v - g)
			return 0;
		--v;
		++l;
		g += ge;
	}
	while (h.good()) {
		if (score == *h - g)
			return 1;
		--h;
		++l;
		g += ge;
	}
	return -1;
}

void traceback(const Sequence &query,
	const Sequence &subject,
	const vector<int> &scores,
	int band,
	int i0,
	int i,
	int j,
	Hsp &l)
{
	Banded_traceback_matrix dp(scores, band, i0);
	//dp.print(i + 1, j + 1);
	l.query_range.end_ = i + 1;
	l.subject_range.end_ = j + 1;
	l.transcript.clear();

	int gap_len, score;

	while ((score = dp(i, j)) > 0) {
		const int match_score = score_matrix(query[i], subject[j]);
		//printf("i=%i j=%i score=%i subject=%c query=%c\n", i, j, dp(i, j), value_traits.alphabet[subject[j]], value_traits.alphabet[query[i]]);

		if (score == match_score + dp(i - 1, j - 1)) {
			if (query[i] == subject[j]) {
				l.transcript.push_back(op_match);
				++l.identities;
				++l.positives;
			}
			else {
				l.transcript.push_back(op_substitution, subject[j]);
				++l.mismatches;
				if (match_score > 0)
					++l.positives;
			}
			--i;
			--j;
			++l.length;
		}
		else {
			const int g = have_gap(dp, i, j, gap_len);
			if (g == -1)
				throw std::runtime_error("Traceback error.");
			++l.gap_openings;
			l.length += gap_len;
			l.gaps += gap_len;
			if (g == 0) {
				i -= gap_len;
				l.transcript.push_back(op_insertion, (unsigned)gap_len);
			}
			else {
				for (; gap_len > 0; gap_len--)
					l.transcript.push_back(op_deletion, subject[j--]);
			}
		}
	}

	l.query_range.begin_ = i + 1;
	l.subject_range.begin_ = j + 1;
	l.transcript.reverse();
	l.transcript.push_terminator();
}

struct Banded_dp_matrix
{

	struct Column_iterator
	{
		inline Column_iterator(const pair<int*, int*> &score, int *hgap) :
			score_(score),
			hgap_(hgap + 1)
		{
		}
		inline int& score()
		{
			return *score_.second;
		}
		inline int diag() const
		{
			return *score_.first;
		}
		inline int hgap_in() const
		{
			return *hgap_;
		}
		inline int& hgap_out()
		{
			return *(hgap_ - 1);
		}
		inline void operator++()
		{
			++score_.first;
			++score_.second;
			++hgap_;
		}
	private:
		pair<int*, int*> score_;
		int *hgap_;
	};

	inline Column_iterator column(int j, int offset)
	{
		return Column_iterator(std::make_pair(&score_[j*band_ + offset], &score_[(j + 1)*band_ + offset]), hgap_.data() + offset);
	}

	const vector<int>& scores() const
	{
		return score_;
	}

	inline Banded_dp_matrix(int band, int cols) :
		band_(band),
		score_(band*(cols + 1)),
		hgap_(band + 1)
	{}

private:

	const int band_;
	vector<int> score_, hgap_;

};

void banded_sw(const Sequence &query, const Sequence &subject, int d_begin, int d_end, int j_begin, int j_end, Hsp &out)
{
	using std::max;
	assert(d_end > d_begin);
	const int slen = (int)subject.length(),
		qlen = (int)query.length();
	d_begin = std::max(d_begin, -(slen - 1));
	d_end = std::min(d_end, qlen);
	const int band = d_end - d_begin,
		gap_open = score_matrix.gap_open() + score_matrix.gap_extend(),
		gap_extend = score_matrix.gap_extend();
	int i0 = d_begin + j_begin, score = 0, max_i, max_j;
	int n=0;
	Banded_dp_matrix mtx(band, slen);
	for (int j = j_begin; j < j_end; ++j, ++i0) {
		const int i1 = std::min(i0 + band, qlen);
		int i = std::max(i0, 0), vgap = 0;
		Banded_dp_matrix::Column_iterator it = mtx.column(j, i - i0);
		for (; i < i1; ++i, ++it) {
			const int match_score = score_matrix(query[i], subject[j]);
			int hgap = it.hgap_in();
			int s = it.diag() + match_score;
			if (s < hgap)
				s = hgap;
			if (s < vgap)
				s = vgap;			
			if (s < 0)
				s = 0;
			const int open = s - gap_open;
			vgap -= gap_extend;
			if (vgap < open)
				vgap = open;
			hgap -= gap_extend;
			if (hgap < open)
				hgap = open;
			it.hgap_out() = hgap;
			it.score() = s;
			if (s > score) {
				score = s;
				max_i = i;
				max_j = j;
			}
		}
	}
	out.score = score;
	traceback(query, subject, mtx.scores(), band, d_begin + j_begin, max_i, max_j, out);
}