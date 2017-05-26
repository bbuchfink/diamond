/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef ALIGN_RANGE_H_
#define ALIGN_RANGE_H_

#include "../basic/statistics.h"
#include "../data/sorted_list.h"
#include "../search/trace_pt_buffer.h"
#include "sse_dist.h"
#include "../dp/dp.h"

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


inline unsigned popcount32(unsigned x)
{
#ifdef _MSC_VER
	return __popcnt(x);
#else
	return __builtin_popcount(x);
#endif
}

inline unsigned popcount64(unsigned long long x)
{
#ifdef _MSC_VER
	return (unsigned)__popcnt64(x);
#else
	return __builtin_popcountll(x);
#endif
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
		vq(TLS::get(vq_ptr)),
		vs(TLS::get(vs_ptr)),
		hits(TLS::get(hits_ptr)),
		stats(stats),
		out(out),
		sid(sid)
	{}
	void run(const sorted_list::const_iterator &q, const sorted_list::const_iterator &s);
	void tiled_search(vector<Finger_print>::const_iterator q,
		vector<Finger_print>::const_iterator q_end,
		vector<Finger_print>::const_iterator s,
		vector<Finger_print>::const_iterator s_end,
		const Range_ref &ref,
		unsigned level);

	static TLS_PTR vector<Finger_print> *vq_ptr, *vs_ptr;
	static TLS_PTR vector<Stage1_hit> *hits_ptr;
	vector<Finger_print> &vq, &vs;
	vector<Stage1_hit> &hits;
	Statistics &stats;
	Trace_pt_buffer::Iterator &out;
	const unsigned sid;
};

void stage2_search(const sorted_list::const_iterator &q,
	const sorted_list::const_iterator &s,
	const vector<Stage1_hit> &hits,
	Statistics &stats,
	Trace_pt_buffer::Iterator &out,
	const unsigned sid);

inline void align_partition(unsigned hp,
		Statistics &stats,
		unsigned sid,
		sorted_list::const_iterator i,
		sorted_list::const_iterator j,
		unsigned thread_id)
{
#ifndef SIMPLE_SEARCH
	if (hp > 0)
		return;
#endif
	Trace_pt_buffer::Iterator* out = new Trace_pt_buffer::Iterator (*Trace_pt_buffer::instance, thread_id);
	Seed_filter seed_filter(stats, *out, sid);
	while(!i.at_end() && !j.at_end()) {
		if(i.key() < j.key()) {
			++i;
		} else if(j.key() < i.key()) {
			++j;
		} else {
			if(i[0] != 0)
				seed_filter.run(j, i);
			++i;
			++j;
		}
	}
	delete out;
}

#endif /* ALIGN_RANGE_H_ */
