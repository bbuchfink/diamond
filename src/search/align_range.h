/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef ALIGN_RANGE_H_
#define ALIGN_RANGE_H_

#include "../basic/statistics.h"
#include "../search/trace_pt_buffer.h"
#include "sse_dist.h"
#include "../dp/dp.h"
#include "../util/intrin.h"
#include "../basic/shape_config.h"

void setup_search_params(pair<size_t, size_t> query_len_bounds, size_t chunk_db_letters);
void setup_search();
void setup_search_cont();

struct Stage1_hit
{
	Stage1_hit(unsigned q_ref, unsigned q_offset, unsigned s_ref, unsigned s_offset) :
		q(q_ref + q_offset),
		s(s_ref + s_offset)
	{}
	bool operator<(const Stage1_hit &rhs) const
	{
		return q < rhs.q;
	}
	struct Query
	{
		unsigned operator()(const Stage1_hit &x) const
		{
			return x.q;
		}
	};
	unsigned q, s;
};

inline int stage2_ungapped(const Letter *query, const Letter *subject, unsigned sid, unsigned &delta, unsigned &len)
{
	return xdrop_ungapped(query, subject, shapes[sid].length_, delta, len);
}

#ifdef __SSE2__

struct Byte_finger_print_48
{
	Byte_finger_print_48(const Letter *q) :
		r1(_mm_loadu_si128((__m128i const*)(q - 16))),
		r2(_mm_loadu_si128((__m128i const*)(q))),
		r3(_mm_loadu_si128((__m128i const*)(q + 16)))
	{}
	static uint64_t match_block(__m128i x, __m128i y)
	{
		return (uint64_t)_mm_movemask_epi8(_mm_cmpeq_epi8(x, y));
	}
	unsigned match(const Byte_finger_print_48 &rhs) const
	{
		return popcount64(match_block(r3, rhs.r3) << 32 | match_block(r1, rhs.r1) << 16 | match_block(r2, rhs.r2));
	}
	__m128i r1, r2, r3;
};

#else

struct Byte_finger_print_48
{
	Byte_finger_print_48()
	{}
	Byte_finger_print_48(const Letter *q)
	{
		//printf("%llx\n", q);
		memcpy(r, q - 16, 48);
	}
	unsigned match(const Byte_finger_print_48 &rhs) const
	{
		unsigned n = 0;
		for (unsigned i = 0; i < 48; ++i)
			if (r[i] == rhs.r[i])
				++n;
		return n;
	}
	Letter r[48];
	//char r[32];
};

#endif

typedef Byte_finger_print_48 Finger_print;

struct Range_ref
{
	Range_ref(vector<Finger_print>::const_iterator q_begin, vector<Finger_print>::const_iterator s_begin) :
		q_begin(q_begin),
		s_begin(s_begin)
	{}
	const vector<Finger_print>::const_iterator q_begin, s_begin;
};

struct Seed_filter
{
	Seed_filter(Statistics &stats, Trace_pt_buffer::Iterator &out, const unsigned sid) :
		stats(stats),
		out(out),
		sid(sid)
	{}
	void run(const Packed_loc *q, size_t nq, const Packed_loc *s, size_t ns);
	void tiled_search(vector<Finger_print>::const_iterator q,
		vector<Finger_print>::const_iterator q_end,
		vector<Finger_print>::const_iterator s,
		vector<Finger_print>::const_iterator s_end,
		const Range_ref &ref,
		unsigned level);

	static thread_local vector<Finger_print> vq, vs;
	static thread_local vector<Stage1_hit> hits;
	Statistics &stats;
	Trace_pt_buffer::Iterator &out;
	const unsigned sid;
};

void stage2_search(const Packed_loc *q,
	const Packed_loc *s,
	const vector<Stage1_hit> &hits,
	Statistics &stats,
	Trace_pt_buffer::Iterator &out,
	const unsigned sid);

#endif /* ALIGN_RANGE_H_ */
