/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#pragma once
#include <stdint.h>
#include <algorithm>
#include "../dp.h"
#include "../../basic/value.h"
#include "../../util/simd/vector.h"
#include "../../stats/hauser_correction.h"

template<typename T, int N>
struct SmallVector
{
	SmallVector() :
		n(0)
	{}
	T& operator[](int i)
	{
		return data[i];
	}
	const T& operator[](int i) const
	{
		return data[i];
	}
	int size() const
	{
		return n;
	}
	void push_back(const T& x)
	{
		data[n++] = x;
	}
	void erase(int i)
	{
		memmove(&data[i], &data[i + 1], (--n - i) * sizeof(T));
	}
private:
	T data[N];
	int n;
};

namespace DISPATCH_ARCH {

template<typename _t>
struct TargetIterator
{

	typedef ::DISPATCH_ARCH::SIMD::Vector<_t> SeqVector;
	enum {
		CHANNELS = SeqVector::CHANNELS
	};

	TargetIterator(std::vector<DpTarget>::const_iterator subject_begin, std::vector<DpTarget>::const_iterator subject_end, int i1, int qlen, int *d_begin) :
		next(0),
		n_targets(int(subject_end - subject_begin)),
		cols(0),
		custom_matrix_16bit(false),
		subject_begin(subject_begin)
	{
		for (; next < std::min((int)CHANNELS, n_targets); ++next) {
			const DpTarget &t = subject_begin[next];
			pos[next] = i1 - (t.d_end - 1);
			const int d0 = d_begin[next];
			//const int d0 = t.d_begin;
			const int j1 = std::min(qlen - 1 - d0, (int)(t.seq.length() - 1)) + 1;
			cols = std::max(cols, j1 - pos[next]);
			target[next] = next;
			active.push_back(next);
			if (t.adjusted_matrix() && (t.matrix->score_max > SCHAR_MAX || t.matrix->score_min < SCHAR_MIN))
				custom_matrix_16bit = true;
		}
	}

	uint64_t live() const {
		uint64_t n = 0;
		for (int i = 0; i < active.size(); ++i) {
			const int channel = active[i];
			if (pos[channel] >= 0)
				n |= 1llu << channel;
		}
		return n;
	}

	char operator[](int channel) const
	{
		if (pos[channel] >= 0) {
			return subject_begin[target[channel]].seq[pos[channel]];
		} else
			return SUPER_HARD_MASK;
	}


	SeqVector get() const
	{
		alignas(32) _t s[CHANNELS];
		std::fill(s, s + CHANNELS, SUPER_HARD_MASK);
		for (int i = 0; i < active.size(); ++i) {
			const int channel = active[i];
			s[channel] = (*this)[channel];
		}
		return SeqVector(s);
	}

	const int8_t** get(const int8_t** target_scores) const {
		static const int8_t blank[32] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		std::fill(target_scores, target_scores + 32, blank);
		for (int i = 0; i < active.size(); ++i) {
			const int channel = active[i];
			const int l = (int)(*this)[channel];
			const DpTarget& dp_target = subject_begin[target[channel]];
			target_scores[channel] = dp_target.adjusted_matrix() ? &dp_target.matrix->scores[32 * l] : &score_matrix.matrix8()[32 * l];
		}
		return target_scores;
	}

	std::vector<const int32_t*> get32() const {
		static const int32_t blank[32] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		std::vector<const int32_t*> target_scores(CHANNELS);
		std::fill(target_scores.begin(), target_scores.end(), blank);
		for (int i = 0; i < active.size(); ++i) {
			const int channel = active[i];
			const int l = (int)(*this)[channel];
			const DpTarget& dp_target = subject_begin[target[channel]];
			target_scores[channel] = dp_target.adjusted_matrix() ? &dp_target.matrix->scores32[32 * l] : &score_matrix.matrix32()[32 * l];
		}
		return target_scores;
	}

	bool init_target(int i, int channel)
	{
		if (next < n_targets) {
			pos[channel] = 0;
			target[channel] = next++;
			return true;
		}
		active.erase(i);
		return false;
	}

	bool inc(int channel)
	{
		++pos[channel];
		if (pos[channel] >= (int)subject_begin[target[channel]].seq.length())
			return false;
		return true;
	}

	uint32_t cbs_mask() const {
		uint32_t r = 0;
		for (uint32_t i = 0; i < (uint32_t)n_targets; ++i)
			if (subject_begin[i].adjusted_matrix())
				r |= 1 << i;
		return r;
	}

	int pos[CHANNELS], target[CHANNELS], next, n_targets, cols;
	bool custom_matrix_16bit;
	SmallVector<int, CHANNELS> active;
	const std::vector<DpTarget>::const_iterator subject_begin;
};

template<typename T, typename It>
struct AsyncTargetBuffer
{

	typedef ::DISPATCH_ARCH::SIMD::Vector<T> SeqVector;
	enum { CHANNELS = SeqVector::CHANNELS };

	AsyncTargetBuffer(const It begin, const It end, std::atomic<BlockId>* const next):
		begin(begin),
		target_count(BlockId(end - begin)),
		next(next),
		custom_matrix_16bit(false)
	{
		BlockId n;
		int i = 0;
		while (i < CHANNELS && (n = (*next)++) < target_count) {
			DpTarget t = begin[n];
			if (t.blank())
				t.target_idx = n;
			pos[i] = 0;
			dp_targets[i] = t;
			active.push_back(i);
			++i;
		}
	}

	int max_len() const {
		int l = 0;
		for (BlockId i = 0; i < target_count; ++i)
			l = std::max(l, (int)DpTarget(begin[i]).seq.length());
		return l;
	}

	char operator[](int channel) const
	{
		if (pos[channel] >= 0) {
			return dp_targets[channel].seq[pos[channel]];
		}
		else
			return SUPER_HARD_MASK;
	}

	SeqVector seq_vector() const
	{
		alignas(32) T s[CHANNELS];
		std::fill(s, s + CHANNELS, SUPER_HARD_MASK);
		for (int i = 0; i < active.size(); ++i) {
			const int channel = active[i];
			s[channel] = (*this)[channel];
		}
		return SeqVector(s);
	}

	const int8_t** get(const int8_t** target_scores) const {
		static const int8_t blank[32] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		std::fill(target_scores, target_scores + 32, blank);
		for (int i = 0; i < active.size(); ++i) {
			const int channel = active[i];
			const int l = (int)(*this)[channel];
			const DpTarget& dp_target = dp_targets[channel];
			target_scores[channel] = dp_target.adjusted_matrix() ? &dp_target.matrix->scores[32 * l] : &score_matrix.matrix8()[32 * l];
		}
		return target_scores;
	}

	std::vector<const int32_t*> get32() const {
		static const int32_t blank[32] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		std::vector<const int32_t*> target_scores(CHANNELS);
		std::fill(target_scores.begin(), target_scores.end(), blank);
		for (int i = 0; i < active.size(); ++i) {
			const int channel = active[i];
			const int l = (int)(*this)[channel];
			const DpTarget& dp_target = dp_targets[channel];
			target_scores[channel] = dp_target.adjusted_matrix() ? &dp_target.matrix->scores32[32 * l] : &score_matrix.matrix32()[32 * l];
		}
		return target_scores;
	}

	bool init_target(int i, int channel)
	{
		const BlockId n = (*next)++;
		if (n >= target_count) {
			active.erase(i);
			return false;
		}
		DpTarget t = begin[n];
		if (t.blank())
			t.target_idx = n;
		pos[channel] = 0;
		dp_targets[channel] = t;
		return true;
	}

	bool inc(int channel)
	{
		++pos[channel];
		if (pos[channel] >= (int)dp_targets[channel].seq.length())
			return false;
		return true;
	}

	uint32_t cbs_mask() {
		uint32_t r = 0;
		custom_matrix_16bit = false;
		for (int i = 0; i < active.size(); ++i) {
			const int channel = active[i];
			if (dp_targets[channel].adjusted_matrix()) {
				r |= 1 << channel;
				if (dp_targets[channel].matrix->score_max > SCHAR_MAX || dp_targets[channel].matrix->score_min < SCHAR_MIN)
					custom_matrix_16bit = true;
			}
		}
		return r;
	}

	int pos[CHANNELS];
	SmallVector<int, CHANNELS> active;
	const It begin;
	const BlockId target_count;
	std::atomic<BlockId>* const next;
	DpTarget dp_targets[CHANNELS];
	bool custom_matrix_16bit;

};

}