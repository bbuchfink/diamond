/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

#include <array>
#include <deque>
#include "basic/value.h"
#include "util/enum.h"
#include "util/intrin.h"

enum struct MaskingAlgo { NONE = 0, TANTAN = 1, SEG = 2, MOTIF = 4 };

DEFINE_ENUM_OPERATORS(MaskingAlgo)

template<>
struct EnumTraits<MaskingAlgo> {
	static const EMap<MaskingAlgo> to_string;
	static const SEMap<MaskingAlgo> from_string;
};

struct MaskingStat {
	MaskingStat() : masked_letters{ 0, 0, 0 } {}
	void add(MaskingAlgo algo, uint64_t n) {
		masked_letters[ctz((uint32_t)algo)] += n;
	}
	uint64_t get(MaskingAlgo algo) const {
		return masked_letters[ctz((uint32_t)algo)];
	}
	MaskingStat& operator+=(const MaskingStat& other) {
		for (size_t i = 0; i < masked_letters.size(); ++i)
			masked_letters[i] += other.masked_letters[i];
		return *this;
	}
	void print(std::ostream& str) const {
		str << "Masked letters: "
			<< "  tantan: " << get(MaskingAlgo::TANTAN)
			<< "  seg: " << get(MaskingAlgo::SEG)
			<< "  motif: " << get(MaskingAlgo::MOTIF) << std::endl;
	}
	std::array<uint64_t, 3> masked_letters;
};

enum class MaskingMode { NONE, TANTAN, BLAST_SEG };

template<>
struct EnumTraits<MaskingMode> {
	static const SEMap<MaskingMode> from_string;
	static const EMap<MaskingMode> to_string;
};

namespace Mask {

struct Ranges : public std::deque<std::pair<Loc, Loc>> {
	void push_back(Loc begin, Loc end) {
		if (empty() || begin > back().second)
			emplace_back(begin, end);
		else
			back().second = end;
	}
	void push_front(Loc loc) {
		if (!empty() && front().first == loc + 1)
			front().first = loc;
		else
			emplace_front(loc, loc + 1);
	}
};

}