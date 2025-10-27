/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <deque>
#include "basic/value.h"
#include "util/enum.h"

enum struct MaskingAlgo { NONE = 0, TANTAN = 1, SEG = 2, MOTIF = 4 };

DEFINE_ENUM_FLAG_OPERATORS(MaskingAlgo)

template<>
struct EnumTraits<MaskingAlgo> {
	static const EMap<MaskingAlgo> to_string;
	static const SEMap<MaskingAlgo> from_string;
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