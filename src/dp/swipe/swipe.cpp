/****
DIAMOND protein aligner
Copyright (C) 2013-2019 Benjamin Buchfink <buchfink@gmail.com>

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
#include "../score_vector.h"
#include "swipe.h"
#include "../../basic/sequence.h"
#include "target_iterator.h"

// #define SW_ENABLE_DEBUG

using std::vector;
using std::pair;

namespace DP { namespace Swipe {
namespace DISPATCH_ARCH {

template<typename _sv>
struct DPMatrix
{
	struct ColumnIterator
	{
		ColumnIterator(_sv* hgap_front, _sv* score_front) :
			hgap_ptr_(hgap_front),
			score_ptr_(score_front)
		{ }
		inline void operator++()
		{
			++hgap_ptr_; ++score_ptr_;
		}
		inline _sv hgap() const
		{
			return *hgap_ptr_;
		}
		inline _sv diag() const
		{
			return *score_ptr_;
		}
		inline void set_hgap(const _sv& x)
		{
			*hgap_ptr_ = x;
		}
		inline void set_score(const _sv& x)
		{
			*score_ptr_ = x;
		}
		_sv *hgap_ptr_, *score_ptr_;
	};
	DPMatrix(int rows)
	{
		hgap_.clear();
		hgap_.resize(rows);
		score_.clear();
		score_.resize(rows + 1);
		std::fill(hgap_.begin(), hgap_.end(), ScoreTraits<_sv>::zero());
		std::fill(score_.begin(), score_.end(), ScoreTraits<_sv>::zero());
	}
	inline ColumnIterator begin()
	{
		return ColumnIterator(&hgap_[0], &score_[0]);
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
	static thread_local std::vector<_sv> hgap_, score_;
};

template<typename _sv> thread_local std::vector<_sv> DPMatrix<_sv>::hgap_;
template<typename _sv> thread_local std::vector<_sv> DPMatrix<_sv>::score_;

#ifdef __SSE2__

template<typename _sv>
vector<int> swipe(const sequence &query, const sequence *subject_begin, const sequence *subject_end)
{
#ifdef SW_ENABLE_DEBUG
	static int v[1024][1024];
#endif

	const int qlen = (int)query.length();
	DPMatrix<_sv> dp(qlen);

	const _sv open_penalty(static_cast<char>(score_matrix.gap_open() + score_matrix.gap_extend())),
		extend_penalty(static_cast<char>(score_matrix.gap_extend())),
		vbias(score_matrix.bias());
	_sv best;
	SwipeProfile<_sv> profile;
	TargetBuffer<_sv::CHANNELS> targets(subject_begin, subject_end);
	vector<int> out(targets.n_targets);

	while (targets.active.size() > 0) {
		typename DPMatrix<_sv>::ColumnIterator it(dp.begin());
		_sv vgap, hgap, last;
		profile.set(targets.seq_vector());
		for (int i = 0; i < qlen; ++i) {
			hgap = it.hgap();
			const _sv next = cell_update<_sv>(it.diag(), profile.get(query[i]), extend_penalty, open_penalty, hgap, vgap, best, vbias);
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

	return out;
}

#endif

vector<int> swipe(const sequence &query, const sequence *subject_begin, const sequence *subject_end)
{
#ifdef __SSE2__
	return swipe<score_vector<uint8_t>>(query, subject_begin, subject_end);
#endif
}

}}}