#pragma once

#include <deque>
#include "../basic/value.h"
#include "../util/enum.h"

enum struct MaskingAlgo { NONE = 0, TANTAN = 1, SEG = 2, MOTIF = 4 };

DEFINE_ENUM_FLAG_OPERATORS(MaskingAlgo)

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