/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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
#include "basic/packed_loc.h"
#include "basic/value.h"

namespace Search {

#ifndef __sparc__
#pragma pack(1)
#endif

struct Hit
{
	using Key = uint32_t;
	using SeedOffset = Loc;

	BlockId query_;
	PackedLoc subject_;
	SeedOffset seed_offset_;
	uint16_t score_;
#ifdef HIT_KEEP_TARGET_ID
	uint32_t target_block_id;
#endif

	Hit() :
		query_(),
		subject_(),
		seed_offset_(),
		score_()
	{ }
	Hit(unsigned query, PackedLoc subject, SeedOffset seed_offset, uint16_t score = 0, uint32_t target_block_id = 0) :
		query_(query),
		subject_(subject),
		seed_offset_(seed_offset),
		score_(score)
#ifdef HIT_KEEP_TARGET_ID
		,target_block_id(target_block_id)
#endif
	{ }
	bool operator==(const Hit& h) const {
		return query_ == h.query_ && subject_ == h.subject_ && seed_offset_ == h.seed_offset_ && score_ == h.score_;
	}
	bool operator<(const Hit& rhs) const
	{
		return query_ < rhs.query_;
	}
	bool blank() const
	{
		return uint64_t(subject_) == 0;
	}
	unsigned operator%(unsigned i) const
	{
		return (query_ / align_mode.query_contexts) % i;
	}
	unsigned operator/(size_t i) const
	{
		return (query_ / align_mode.query_contexts) / (unsigned)i;
	}
	unsigned frame() const {
		return query_ % align_mode.query_contexts;
	}
	int64_t global_diagonal() const
	{
		return (int64_t)subject_ - (int64_t)seed_offset_;
	}
	template<unsigned d>
	static unsigned query_id(const Hit& x)
	{
		return x.query_ / d;
	}
	template<unsigned d>
	struct Query_id
	{
		unsigned operator()(const Hit& x) const
		{
			return query_id<d>(x);
		}
	};
	struct Query {
		unsigned operator()(const Hit& h) const {
			return h.query_;
		}
	};
	struct SourceQuery {
		int32_t operator()(const Hit& h) const {
			return h.query_ / contexts;
		}
		const int32_t contexts;
	};
	struct Subject {
		uint64_t operator()(const Hit& h) const {
			return h.subject_;
		}
	};
	struct CmpSubject {
		bool operator()(const Hit& lhs, const Hit& rhs) const
		{
			return lhs.subject_ < rhs.subject_
				|| (lhs.subject_ == rhs.subject_ && (lhs.query_ < rhs.query_ || (lhs.query_ == rhs.query_ && lhs.seed_offset_ < rhs.seed_offset_)));
		}
	};
	struct CmpQueryTarget {
		bool operator()(const Hit& x, const Hit& y) const {
			return x.query_ < y.query_ || (x.query_ == y.query_ && x.subject_ < y.subject_);
		}
	};
	struct CmpTargetOffset {
		bool operator()(const Hit& x, size_t s) const {
			return (uint64_t)x.subject_ < s;
		}
	};
	static bool cmp_normalized_subject(const Hit &lhs, const Hit &rhs)
	{
		const uint64_t x = (uint64_t)lhs.subject_ + (uint64_t)rhs.seed_offset_, y = (uint64_t)rhs.subject_ + (uint64_t)lhs.seed_offset_;
		return x < y || (x == y && lhs.seed_offset_ < rhs.seed_offset_);
	}
	static bool cmp_frame(const Hit &x, const Hit &y) {
		return x.frame() < y.frame();
	}
	friend std::ostream& operator<<(std::ostream &s, const Hit &me)
	{
		s << me.query_ << '\t' << uint64_t(me.subject_) << '\t' << me.seed_offset_ << '\t' << me.score_ << '\n';
		return s;
	}
} PACKED_ATTRIBUTE;

#pragma pack()

}