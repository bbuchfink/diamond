/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <array>
#include <deque>
#include "basic/value.h"
#include "util/enum.h"
#include "util/intrin.h"
#include "util/log_stream.h"

enum struct MaskingAlgo { NONE = 0, TANTAN = 1, SEG = 2, MOTIF = 4 };

DEFINE_ENUM_FLAG_OPERATORS(MaskingAlgo)

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
	void print(MessageStream& str) const {
		str << "Masked letters: "
			<< "  tantan: " << get(MaskingAlgo::TANTAN)
			<< "  seg: " << get(MaskingAlgo::SEG)
			<< "  motif: " << get(MaskingAlgo::MOTIF) << std::endl;
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