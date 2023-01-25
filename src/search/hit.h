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
#include <stdint.h>
#include "../basic/packed_loc.h"
#include "../basic/value.h"
#include "../util/io/input_file.h"
#include "../util/io/serialize.h"

namespace Search {

#pragma pack(1)

struct Hit
{
	using Key = uint32_t;
	using SeedOffset = Loc;

	uint32_t query_;
	PackedLoc subject_;
	SeedOffset seed_offset_;
	uint16_t score_;
#ifdef HIT_KEEP_TARGET_ID
	uint32_t target_block_id;
#endif

	Hit() :
		query_(),
		subject_(),
		seed_offset_()
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
	template<unsigned _d>
	static unsigned query_id(const Hit& x)
	{
		return x.query_ / _d;
	}
	template<unsigned _d>
	struct Query_id
	{
		unsigned operator()(const Hit& x) const
		{
			return query_id<_d>(x);
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
				|| (lhs.subject_ == rhs.subject_ && lhs.seed_offset_ < rhs.seed_offset_);
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

template<> struct SerializerTraits<Search::Hit> {
	using Key = BlockId;
	SerializerTraits(bool long_subject_offsets, int32_t query_contexts):
		long_subject_offsets(long_subject_offsets),
		key{ query_contexts }
	{}
	const bool long_subject_offsets;
	const struct {
		Key operator()(const Search::Hit& hit) const {
			return hit.query_ / query_contexts;
		}
		const int32_t query_contexts;
	} key;
	static Search::Hit make_sentry(uint32_t query, Loc seed_offset) {
		return { query, 0, seed_offset,0 };
	}
	static bool is_sentry(const Search::Hit& hit) {
		return hit.score_ == 0;
	}
};

template<> struct TypeSerializer<Search::Hit> {

	TypeSerializer(TextBuffer& buf, const SerializerTraits<Search::Hit>& traits):
		traits(traits),
		buf_(&buf)
	{}

	TypeSerializer& operator<<(const Search::Hit& hit) {
		if (SerializerTraits<Search::Hit>::is_sentry(hit)) {
			buf_->write((uint16_t)0);
			buf_->write_varint(hit.query_);
			buf_->write_varint(hit.seed_offset_);
			return *this;
		}
		buf_->write((uint16_t)hit.score_);
		if (traits.long_subject_offsets)
			buf_->write_raw((const char*)&hit.subject_, 5);
		else
			buf_->write(hit.subject_.low);
#ifdef HIT_KEEP_TARGET_ID
		buf_->write(hit.target_block_id);
#endif
		return *this;
	}

	const SerializerTraits<Search::Hit> traits;

private:

	TextBuffer* buf_;

};

template<> struct TypeDeserializer<Search::Hit> {

	TypeDeserializer(InputFile& f, const SerializerTraits<Search::Hit>& traits):
		f_(&f),
		traits_(traits)
	{
	}

	template<typename It>
	TypeDeserializer<Search::Hit>& operator>>(It& it) {
		uint16_t x;
		f_->read(x);

		for (;;) {
			uint32_t query_id, seed_offset;
			f_->varint = true;
			(*f_) >> query_id >> seed_offset;
			PackedLoc subject_loc;
			uint32_t x;
			f_->varint = false;
			for (;;) {
				uint16_t score;
				try {
					f_->read(score);
				}
				catch (EndOfStream&) {
					return *this;
				}
				if (score == 0)
					break;
				if (traits_.long_subject_offsets)
					f_->read(subject_loc);
				else {
					f_->read(x);
					subject_loc = x;
				}
#ifdef HIT_KEEP_TARGET_ID
				uint32_t target_block_id;
				f_->read(target_block_id);
				*it = { query_id, subject_loc, (Loc)seed_offset, score, target_block_id };
#else
				*it = { query_id, subject_loc, (Loc)seed_offset, score };
#endif
			}
		}
	}

private:

	InputFile* f_;
	const SerializerTraits<Search::Hit> traits_;

};
