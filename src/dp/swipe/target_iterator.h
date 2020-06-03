/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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
#include "../basic/value.h"
#include "../../util/simd/vector.h"

namespace DISPATCH_ARCH {

template<typename _t>
struct TargetIterator
{

	typedef ::DISPATCH_ARCH::SIMD::Vector<_t> SeqVector;
	enum {
		CHANNELS = SeqVector::CHANNELS
	};

	TargetIterator(vector<DpTarget>::const_iterator subject_begin, vector<DpTarget>::const_iterator subject_end, int i1, int qlen, int *d_begin) :
		next(0),
		n_targets(int(subject_end - subject_begin)),
		cols(0),
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

#ifdef __SSSE3__
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
#else
	uint64_t get() const
	{
		uint64_t dst = 0;
		for (int i = 0; i < active.size(); ++i) {
			const int channel = active[i];
			dst |= uint64_t((*this)[channel]) << (8 * channel);
		}
		return dst;
	}
#endif

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

	int pos[CHANNELS], target[CHANNELS], next, n_targets, cols;
	Static_vector<int, CHANNELS> active;
	const vector<DpTarget>::const_iterator subject_begin;
};

template<typename _t>
struct TargetBuffer
{

	typedef ::DISPATCH_ARCH::SIMD::Vector<_t> SeqVector;
	enum { CHANNELS = SeqVector::CHANNELS };

	TargetBuffer(const sequence *subject_begin, const sequence *subject_end) :
		next(0),
		n_targets(int(subject_end - subject_begin)),
		subject_begin(subject_begin)
	{
		for (; next < std::min((int)CHANNELS, n_targets); ++next) {
			pos[next] = 0;
			target[next] = next;
			active.push_back(next);
		}
	}

	char operator[](int channel) const
	{
		if (pos[channel] >= 0) {
			return subject_begin[target[channel]][pos[channel]];
		}
		else
			return value_traits.mask_char;
	}

#ifdef __SSSE3__
	SeqVector seq_vector() const
	{
		alignas(32) _t s[CHANNELS];
		std::fill(s, s + CHANNELS, SUPER_HARD_MASK);
		for (int i = 0; i < active.size(); ++i) {
			const int channel = active[i];
			s[channel] = (*this)[channel];
		}
		return SeqVector(s);
	}
#else
	uint64_t seq_vector()
	{
		uint64_t dst = 0;
		for (int i = 0; i < active.size(); ++i) {
			const int channel = active[i];
			dst |= uint64_t((*this)[channel]) << (8 * channel);
		}
		return dst;
	}
#endif

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
		if (pos[channel] >= (int)subject_begin[target[channel]].length())
			return false;
		return true;
	}

	int pos[CHANNELS], target[CHANNELS], next, n_targets, cols;
	Static_vector<int, CHANNELS> active;
	const sequence *subject_begin;
};

}