/****
DIAMOND protein aligner
Copyright (C) 2021 Max Planck Society for the Advancement of Science e.V.

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
#include <ostream>
#include "../../basic/match.h"

FORCE_INLINE uint8_t cmp_mask(int x, int y) {
	return x == y;
}

FORCE_INLINE int blend(int v, int w, int mask) {
	return mask ? w : v;
}

template<typename Sv>
struct DummyIdMask {
	DummyIdMask(const Letter q, const Sv& t)
	{}
};

template<typename Sv>
struct VectorIdMask {
	VectorIdMask(const Letter q, const Sv& t) :
		mask(blend(Sv(0), Sv(1), Sv(typename ::DISPATCH_ARCH::ScoreTraits<Sv>::Score(q)) == t))
	{}
	const Sv mask;
};

template<typename Sv>
struct ForwardCell : public Sv {
	ForwardCell() :
		Sv(),
		ident(),
		len()
	{}
	ForwardCell(const Sv& v) :
		Sv(v),
		ident(),
		len()
	{}
	Sv ident, len;
};

template<>
struct ForwardCell<int32_t> {
	int32_t v, ident, len;
	operator int32_t() const {
		return v;
	}
	ForwardCell(const int32_t v) :
		v(v),
		ident(0),
		len(0)
	{}
	ForwardCell() :
		v(0),
		ident(0),
		len(0)
	{}
	ForwardCell& operator-=(int32_t x) {
		v -= x;
		return *this;
	}
	ForwardCell& operator+=(int32_t x) {
		v += x;
		return *this;
	}
	void max(const ForwardCell& x) {
		v = std::max(v, x.v);
	}
	struct Stats {
		int ident, len;
		friend std::ostream& operator<<(std::ostream& s, const Stats& c) {
			s << " ident=" << c.ident << " len=" << c.len;
			return s;
		}
	};
};

template<typename Sv>
FORCE_INLINE void set_channel(ForwardCell<Sv>& v, const int i, const typename ::DISPATCH_ARCH::ScoreTraits<Sv>::Score x) {
	set_channel((Sv&)v, i, x);
	set_channel(v.ident, i, x);
	set_channel(v.len, i, x);
}

template<typename Sv>
struct BackwardCell : public Sv {
	BackwardCell() :
		Sv(),
		mismatch(),
		gapopen()
	{}
	BackwardCell(const Sv& v) :
		Sv(v),
		mismatch(),
		gapopen()
	{}
	Sv mismatch, gapopen;
};

template<>
struct BackwardCell<int32_t> {
	int32_t v, mismatch, gapopen;
	operator int32_t() const {
		return v;
	}
	BackwardCell(const int32_t v) :
		v(v),
		mismatch(0),
		gapopen(0)
	{}
	BackwardCell() :
		v(0),
		mismatch(0),
		gapopen(0)
	{}
	BackwardCell& operator-=(int32_t x) {
		v -= x;
		return *this;
	}
	BackwardCell& operator+=(int32_t x) {
		v += x;
		return *this;
	}
	void max(const BackwardCell& x) {
		v = std::max(v, x.v);
	}
	struct Stats {
		int mismatch, gap_open;
		friend std::ostream& operator<<(std::ostream& s, const Stats& c) {
			s << " mismatch=" << c.mismatch << " gapopen=" << c.gap_open;
			return s;
		}
	};
};

template<typename Sv>
FORCE_INLINE void set_channel(BackwardCell<Sv>& v, const int i, const typename ::DISPATCH_ARCH::ScoreTraits<Sv>::Score x) {
	set_channel((Sv&)v, i, x);
	set_channel(v.mismatch, i, x);
	set_channel(v.gapopen, i, x);
}

struct Void {
	friend std::ostream& operator<<(std::ostream& s, const Void&) {
		return s;
	}
};

template<typename Sv>
FORCE_INLINE Void extract_stats(const Sv&, int) {
	return Void();
}

template<typename Sv>
FORCE_INLINE ForwardCell<int32_t>::Stats extract_stats(const ForwardCell<Sv>& v, int channel) {
	const auto s = ::DISPATCH_ARCH::ScoreTraits<Sv>::int_score;
	return { s(extract_channel(v.ident, channel)), s(extract_channel(v.len, channel)) };
}

template<typename Sv>
FORCE_INLINE BackwardCell<int32_t>::Stats extract_stats(const BackwardCell<Sv>& v, int channel) {
	const auto s = ::DISPATCH_ARCH::ScoreTraits<Sv>::int_score;
	return { s(extract_channel(v.mismatch, channel)), s(extract_channel(v.gapopen, channel)) };
}

template<typename Sv>
FORCE_INLINE bool overflow_stats(Void) {
	return false;
}

template<typename Sv>
FORCE_INLINE bool overflow_stats(const ForwardCell<int32_t>::Stats& stats) {
	constexpr auto m = DISPATCH_ARCH::ScoreTraits<Sv>::max_int_score();
	return stats.ident == m || stats.len == m;
}

template<typename Sv>
FORCE_INLINE bool overflow_stats(const BackwardCell<int32_t>::Stats& stats) {
	constexpr auto m = DISPATCH_ARCH::ScoreTraits<Sv>::max_int_score();
	return stats.gap_open == m || stats.mismatch == m;
}

FORCE_INLINE void assign_stats(Hsp& hsp, Void) {}

FORCE_INLINE void assign_stats(Hsp& hsp, const ForwardCell<int32_t>::Stats& v) {
	hsp.identities = v.ident;
	hsp.length = v.len;
}

FORCE_INLINE void assign_stats(Hsp& hsp, const BackwardCell<int32_t>::Stats& v) {
	hsp.gap_openings = v.gap_open;
	hsp.mismatches = v.mismatch;
	hsp.gaps = hsp.length - hsp.identities - hsp.mismatches;
}

template<typename Sv>
FORCE_INLINE void update_stats(const Sv&, const Sv&, const Sv&, const DummyIdMask<Sv>&) {
}

template<typename Sv>
FORCE_INLINE void update_stats(ForwardCell<Sv>& current_cell, ForwardCell<Sv>& horizontal_gap, ForwardCell<Sv>& vertical_gap, const VectorIdMask<Sv>& id_mask) {
	const Sv one = Sv(1);
	current_cell.ident += id_mask.mask;
	current_cell.len += one;
	horizontal_gap.len += one;
	vertical_gap.len += one;
}

template<typename Sv>
FORCE_INLINE void update_stats(BackwardCell<Sv>& current_cell, BackwardCell<Sv>& horizontal_gap, BackwardCell<Sv>& vertical_gap, const VectorIdMask<Sv>& id_mask) {
	current_cell.mismatch += Sv(1) - id_mask.mask;
}

template<typename Sv>
FORCE_INLINE void update_open(const Sv&, const Sv&) {
}

template<typename Sv>
FORCE_INLINE void update_open(ForwardCell<Sv>& open, ForwardCell<Sv>& current) {
	const Sv zero = Sv(), zero_mask = current == zero;
	current.ident = blend(current.ident, zero, zero_mask);
	current.len = blend(current.len, zero, zero_mask);
}

template<typename Sv>
FORCE_INLINE void update_open(BackwardCell<Sv>& open, BackwardCell<Sv>& current) {
	open.gapopen += Sv(1);
	const Sv zero = Sv(), zero_mask = current == zero;
	current.mismatch = blend(current.mismatch, zero, zero_mask);
	current.gapopen = blend(current.gapopen, zero, zero_mask);
}

template<typename Sv>
FORCE_INLINE void set_max(ForwardCell<Sv>& v, const ForwardCell<Sv>& x) {
	v.max(x);
	const Sv mask = v == x;
	v.ident = blend(v.ident, x.ident, mask);
	v.len = blend(v.len, x.len, mask);
}

template<typename Sv>
FORCE_INLINE void set_max(BackwardCell<Sv>& v, const BackwardCell<Sv>& x) {
	v.max(x);
	const Sv mask = v == x;
	v.mismatch = blend(v.mismatch, x.mismatch, mask);
	v.gapopen = blend(v.gapopen, x.gapopen, mask);
}

FORCE_INLINE void saturate(ForwardCell<int32_t>& c) {
	c.v = std::max(c.v, 0);
}

FORCE_INLINE void saturate(BackwardCell<int32_t>& c) {
	c.v = std::max(c.v, 0);
}