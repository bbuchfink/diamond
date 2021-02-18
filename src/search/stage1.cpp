/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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

#include <vector>
#include "search.h"
#include "sse_dist.h"
#include "seed_complexity.h"
#include "finger_print.h"
#include "../data/reference.h"
#include "../data/queries.h"
#include "../util/memory/alignment.h"

using std::vector;

namespace Search {
namespace DISPATCH_ARCH {

const unsigned tile_size[] = { 1024, 128 };

constexpr ptrdiff_t INNER_LOOP_QUERIES = 6;
typedef vector<Finger_print, Util::Memory::AlignmentAllocator<Finger_print, 16>> Container;
typedef Container::const_iterator Ptr;

struct Range_ref
{
	Range_ref(Ptr q_begin, Ptr s_begin) :
		q_begin(q_begin),
		s_begin(s_begin)
	{}
	const Ptr q_begin, s_begin;
};

#define FAST_COMPARE2(q, s, stats, q_ref, s_ref, q_offset, s_offset, hits) if (q.match(s) >= config.min_identities) stats.inc(Statistics::TENTATIVE_MATCHES1)
#define FAST_COMPARE(q, s, stats, q_ref, s_ref, q_offset, s_offset, hits) if (q.match(s) >= config.min_identities) hits.push_back(Stage1_hit(q_ref, q_offset, s_ref, s_offset))

void query_register_search(Ptr q,
	Ptr s,
	Ptr s_end,
	const Range_ref &ref,
	vector<Stage1_hit> &hits,
	Statistics &stats)
{
	const unsigned q_ref = unsigned(q - ref.q_begin);
	unsigned s_ref = unsigned(s - ref.s_begin);
	alignas(16) Finger_print q1 = *(q++), q2 = *(q++), q3 = *(q++), q4 = *(q++), q5 = *(q++), q6 = *q;
	const Ptr end2 = s_end - (s_end - s) % 4;
	for (; s < end2; ) {
		alignas(16) Finger_print s1 = *(s++), s2 = *(s++), s3 = *(s++), s4 = *(s++);
		stats.inc(Statistics::SEED_HITS, 6 * 4);
		FAST_COMPARE(q1, s1, stats, q_ref, s_ref, 0, 0, hits);
		FAST_COMPARE(q2, s1, stats, q_ref, s_ref, 1, 0, hits);
		FAST_COMPARE(q3, s1, stats, q_ref, s_ref, 2, 0, hits);
		FAST_COMPARE(q4, s1, stats, q_ref, s_ref, 3, 0, hits);
		FAST_COMPARE(q5, s1, stats, q_ref, s_ref, 4, 0, hits);
		FAST_COMPARE(q6, s1, stats, q_ref, s_ref, 5, 0, hits);
		FAST_COMPARE(q1, s2, stats, q_ref, s_ref, 0, 1, hits);
		FAST_COMPARE(q2, s2, stats, q_ref, s_ref, 1, 1, hits);
		FAST_COMPARE(q3, s2, stats, q_ref, s_ref, 2, 1, hits);
		FAST_COMPARE(q4, s2, stats, q_ref, s_ref, 3, 1, hits);
		FAST_COMPARE(q5, s2, stats, q_ref, s_ref, 4, 1, hits);
		FAST_COMPARE(q6, s2, stats, q_ref, s_ref, 5, 1, hits);
		FAST_COMPARE(q1, s3, stats, q_ref, s_ref, 0, 2, hits);
		FAST_COMPARE(q2, s3, stats, q_ref, s_ref, 1, 2, hits);
		FAST_COMPARE(q3, s3, stats, q_ref, s_ref, 2, 2, hits);
		FAST_COMPARE(q4, s3, stats, q_ref, s_ref, 3, 2, hits);
		FAST_COMPARE(q5, s3, stats, q_ref, s_ref, 4, 2, hits);
		FAST_COMPARE(q6, s3, stats, q_ref, s_ref, 5, 2, hits);
		FAST_COMPARE(q1, s4, stats, q_ref, s_ref, 0, 3, hits);
		FAST_COMPARE(q2, s4, stats, q_ref, s_ref, 1, 3, hits);
		FAST_COMPARE(q3, s4, stats, q_ref, s_ref, 2, 3, hits);
		FAST_COMPARE(q4, s4, stats, q_ref, s_ref, 3, 3, hits);
		FAST_COMPARE(q5, s4, stats, q_ref, s_ref, 4, 3, hits);
		FAST_COMPARE(q6, s4, stats, q_ref, s_ref, 5, 3, hits);
		s_ref += 4;
	}
	for (; s < s_end; ++s) {
		stats.inc(Statistics::SEED_HITS, 6);
		FAST_COMPARE(q1, *s, stats, q_ref, s_ref, 0, 0, hits);
		FAST_COMPARE(q2, *s, stats, q_ref, s_ref, 1, 0, hits);
		FAST_COMPARE(q3, *s, stats, q_ref, s_ref, 2, 0, hits);
		FAST_COMPARE(q4, *s, stats, q_ref, s_ref, 3, 0, hits);
		FAST_COMPARE(q5, *s, stats, q_ref, s_ref, 4, 0, hits);
		FAST_COMPARE(q6, *s, stats, q_ref, s_ref, 5, 0, hits);
		++s_ref;
	}
}

void inner_search(Ptr q,
	Ptr q_end,
	Ptr s,
	Ptr s_end,
	const Range_ref &ref,
	vector<Stage1_hit> &hits,
	Statistics &stats)
{
	unsigned q_ref = unsigned(q - ref.q_begin);
	for (; q < q_end; ++q) {
		unsigned s_ref = unsigned(s - ref.s_begin);
		for (Ptr s2 = s; s2 < s_end; ++s2) {
			stats.inc(Statistics::SEED_HITS);
			FAST_COMPARE((*q), *s2, stats, q_ref, s_ref, 0, 0, hits);
			++s_ref;
		}
		++q_ref;
	}
}

void tiled_search(Ptr q,
	Ptr q_end,
	Ptr s,
	Ptr s_end,
	const Range_ref &ref,
	unsigned level,
	std::vector<Stage1_hit> &hits,
	Statistics &stats)
{
	switch (level) {
	case 0:
	case 1:
		for (; q < q_end; q += std::min(q_end - q, (ptrdiff_t)tile_size[level]))
			for (Ptr s2 = s; s2 < s_end; s2 += std::min(s_end - s2, (ptrdiff_t)tile_size[level]))
				tiled_search(q, q + std::min(q_end - q, (ptrdiff_t)tile_size[level]), s2, s2 + std::min(s_end - s2, (ptrdiff_t)tile_size[level]), ref, level + 1, hits, stats);
		break;
	case 2:
		for (; q < q_end; q += std::min(q_end - q, (ptrdiff_t)INNER_LOOP_QUERIES))
			if (q_end - q < INNER_LOOP_QUERIES)
				inner_search(q, q_end, s, s_end, ref, hits, stats);
			else
				query_register_search(q, s, s_end, ref, hits, stats);
	}
}

static void load_fps(const Packed_loc *p, size_t n, Container &v, const SequenceSet &seqs)
{
	v.clear();
	v.reserve(n);
	const Packed_loc *end = p + n;
	for (; p < end; ++p)
		v.emplace_back(seqs.data(*p));
}

}}