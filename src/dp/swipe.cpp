/****
Copyright (c) 2017, Benjamin Buchfink
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

#include <vector>
#include "dp.h"
#include "score_vector.h"

using std::vector;
using std::pair;

#ifdef __SSE2__

template<typename _score>
struct Swipe_profile
{
	inline void set(const __m128i &seq)
	{
		assert(sizeof(data_) / sizeof(score_vector<_score>) >= value_traits.alphabet_size);
		for (unsigned j = 0; j < value_traits.alphabet_size; ++j)
			data_[j] = score_vector<_score>(j, seq);
	}
	inline const score_vector<_score>& get(Letter i) const
	{
		return data_[(int)i];
	}
	score_vector<_score> data_[25];
};

template<typename _score>
struct Swipe_matrix
{
	typedef score_vector<_score> sv;
	struct Column_iterator
	{
		Column_iterator(sv* hgap_front, sv* score_front) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front)
		{ }
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_;
		}
		inline sv hgap() const
		{
			return *hgap_ptr_;
		}
		inline sv diag() const
		{
			return *score_ptr_;
		}
		inline void set_hgap(const sv& x)
		{
			*hgap_ptr_ = x;
		}
		inline void set_score(const sv& x)
		{
			*score_ptr_ = x;
		}
		sv *hgap_ptr_, *score_ptr_;
	};
	Swipe_matrix(int rows):
		hgap_(TLS::get(hgap_ptr)),
		score_(TLS::get(score_ptr))
	{
		hgap_.clear();
		hgap_.resize(rows);
		score_.clear();
		score_.resize(rows + 1);
		memset(hgap_.data(), 0, rows * sizeof(sv));
		memset(score_.data(), 0, (rows + 1) * sizeof(sv));
	}
	inline Column_iterator begin()
	{
		return Column_iterator(&hgap_[0], &score_[0]);
	}
	void set_zero(int c)
	{
		const int l = (int)hgap_.size();
		for (int i = 0; i < l; ++i) {
			hgap_[i].set(c, 0);
			score_[i].set(c, 0);
		}
		score_[l].set(c, 0);
	}
private:
	vector<sv> &hgap_, &score_;
	static TLS_PTR vector<sv> *hgap_ptr, *score_ptr;
};

template<typename _score> TLS_PTR vector<score_vector<_score> >* Swipe_matrix<_score>::hgap_ptr;
template<typename _score> TLS_PTR vector<score_vector<_score> >* Swipe_matrix<_score>::score_ptr;

template<typename _score>
inline score_vector<_score> cell_update(const score_vector<_score> &diagonal_cell,
	const score_vector<_score> &scores,
	const score_vector<_score> &gap_extension,
	const score_vector<_score> &gap_open,
	score_vector<_score> &horizontal_gap,
	score_vector<_score> &vertical_gap,
	score_vector<_score> &best,
	const score_vector<_score> &vbias)
{
	score_vector<_score> current_cell = diagonal_cell + scores;
	current_cell -= vbias;
	current_cell.max(vertical_gap).max(horizontal_gap);
	best.max(current_cell);
	vertical_gap -= gap_extension;
	horizontal_gap -= gap_extension;
	const score_vector<_score> open = current_cell - gap_open;
	vertical_gap.max(open);
	horizontal_gap.max(open);
	return current_cell;
}

template<int _n>
struct Target_iterator
{
	Target_iterator(vector<sequence>::const_iterator subject_begin, vector<sequence>::const_iterator subject_end):
		next(0),
		n_targets(int(subject_end-subject_begin)),
		subject_begin(subject_begin)
	{
		for (; next < std::min(_n, n_targets); ++next) {
			pos[next] = 0;
			target[next] = next;
			active.push_back(next);
		}
	}
	char operator[](int i) const
	{
		return subject_begin[target[i]][pos[i]];
	}
	__m128i get() const
	{
		char s[16];
		for (int i = 0; i < active.size(); ++i) {
			const int j = active[i];
			s[j] = (*this)[j];
		}
		return _mm_loadu_si128((const __m128i*)s);
	}
	bool init_target(int i, int j)
	{
		if (next < n_targets) {
			pos[j] = 0;
			target[j] = next++;
			return true;
		}
		active.erase(i);
		return false;
	}
	bool inc(int i)
	{
		++pos[i];
		if (pos[i] >= subject_begin[target[i]].length())
			return false;
		return true;
	}
	int pos[_n], target[_n], next, n_targets;
	Static_vector<int, _n> active;
	const vector<sequence>::const_iterator subject_begin;
};

template<typename _score>
void swipe(const sequence &query, vector<sequence>::const_iterator subject_begin, vector<sequence>::const_iterator subject_end, vector<int>::iterator out)
{
#ifdef SW_ENABLE_DEBUG
	int v[1024][1024];
#endif

	typedef score_vector<_score> sv;

	const int qlen = (int)query.length();
	Swipe_matrix<_score> dp(qlen);

	const sv open_penalty(static_cast<char>(score_matrix.gap_open() + score_matrix.gap_extend())),
		extend_penalty(static_cast<char>(score_matrix.gap_extend())),
		vbias(score_matrix.bias());
	sv best;
	Swipe_profile<_score> profile;
	Target_iterator<score_traits<_score>::channels> targets(subject_begin, subject_end);

	while (targets.active.size() > 0) {
		typename Swipe_matrix<_score>::Column_iterator it(dp.begin());
		sv vgap, hgap, last;
		profile.set(targets.get());
		for (int i = 0; i < qlen; ++i) {
			hgap = it.hgap();
			const sv next = cell_update<_score>(it.diag(), profile.get(query[i]), extend_penalty, open_penalty, hgap, vgap, best, vbias);
			it.set_hgap(hgap);
			it.set_score(last);
			last = next;
#ifdef SW_ENABLE_DEBUG
			v[j][it.row_pos_] = next[0];
#endif
			++it;
		}
		it.set_score(last);
		
		for (int i = 0; i < targets.active.size(); ++i) {
			int j = targets.active[i];
			if (!targets.inc(j)) {
				out[targets.target[j]] = best[j];
				if (targets.init_target(i, j)) {
					dp.set_zero(j);
					best.set(j, 0);
				}
			}
		}
	}

#ifdef SW_ENABLE_DEBUG
	for (unsigned j = 0; j<qlen; ++j) {
		for (unsigned i = 0; i<subjects[0].length(); ++i)
			printf("%4i", v[i][j]);
		printf("\n");
	}
	printf("\n");
#endif
}

void swipe(const sequence &query, vector<sequence>::const_iterator subject_begin, vector<sequence>::const_iterator subject_end, vector<int>::iterator out)
{
	swipe<uint8_t>(query, subject_begin, subject_end, out);
}

#endif