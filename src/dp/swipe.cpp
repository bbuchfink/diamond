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

#include <vector>
#include "dp.h"
#include "score_vector.h"

// #define SW_ENABLE_DEBUG

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
		if (pos[i] >= (int)subject_begin[target[i]].length())
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
	static int v[1024][1024];
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
			v[targets.pos[0]][i] = next[0];
#endif
			++it;
		}
		it.set_score(last);
		
		for (int i = 0; i < targets.active.size();) {
			int j = targets.active[i];
			if (!targets.inc(j)) {
				out[targets.target[j]] = best[j];
				if (targets.init_target(i, j)) {
					dp.set_zero(j);
					best.set(j, 0);
				}
				else
					continue;
			}
			++i;
		}
	}

#ifdef SW_ENABLE_DEBUG
	for (unsigned j = 0; j < qlen; ++j) {
		for (unsigned i = 0; i < subject_begin[0].length(); ++i)
			printf("%4i", v[i][j]);
		printf("\n");
	}
	printf("\n");
#endif
}

#endif

void swipe(const sequence &query, vector<sequence>::const_iterator subject_begin, vector<sequence>::const_iterator subject_end, vector<int>::iterator out)
{
#ifdef __SSE2__
	swipe<uint8_t>(query, subject_begin, subject_end, out);
#endif
}
