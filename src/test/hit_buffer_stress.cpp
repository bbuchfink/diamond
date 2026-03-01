/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

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

// Stress test for HitBuffer::load() / retrieve() pipeline.
// Covers both the in-memory (trace_pt_membuf) and disk-backed code paths.

#include <vector>
#include <thread>
#include <atomic>
#include <iostream>
#include <limits>
#include <stdexcept>
#include "basic/config.h"
#include "search/hit_buffer.h"
#include "basic/config.h"
#include "util/parallel/simple_thread_pool.h"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Per-hit fingerprint: order-independent accumulation is safe because we use
// plain addition as the outer reduce operator.
static uint64_t hit_fp(uint32_t query, uint64_t subject, uint32_t seed_offset, uint16_t score)
{
	return (uint64_t)query       * 2654435769ULL
		 + (uint64_t)subject     * 2246822519ULL
		 + (uint64_t)seed_offset * 3266489917ULL
		 + (uint64_t)score       *  668265261ULL;
}

// ---------------------------------------------------------------------------
// Single-mode test
// ---------------------------------------------------------------------------

struct HitBufferTestResult {
	bool   passed;
	size_t expected_hits;
	size_t actual_hits;
	uint64_t expected_checksum;
	uint64_t actual_checksum;
};

static HitBufferTestResult run_single_mode(bool membuf_mode)
{
	// ---- Parameters -------------------------------------------------------
	static constexpr int      BIN_COUNT       = 32;
	static constexpr int      QUERIES_PER_BIN = 100000;
	static constexpr int      QUERY_COUNT     = BIN_COUNT * QUERIES_PER_BIN;
	static constexpr int      HITS_PER_QUERY  = 111;
	static constexpr uint64_t TARGET_LEN      = 50000; // subject locs drawn from [1, TARGET_LEN]
	static constexpr int      THREAD_COUNT    = BIN_COUNT; // one writer thread per bin
	static constexpr int      QUERY_CONTEXTS  = 1;

	// max_query / max_target bounds used by load_bin validation
	const uint32_t max_query  = (uint32_t)QUERY_COUNT;
	const uint64_t max_target = TARGET_LEN + 1; // subject values top out at TARGET_LEN

	// ---- Set global config flags for this run -----------------------------
	const bool saved_membuf = config.trace_pt_membuf;
	const bool saved_swipe  = config.swipe_all;
	config.trace_pt_membuf  = membuf_mode;
	config.swipe_all        = false;

	// ---- Build key_partition: 4 equal bins --------------------------------
	std::vector<uint32_t> key_partition;
	for (int i = 1; i <= BIN_COUNT; ++i)
		key_partition.push_back(uint32_t(i * QUERIES_PER_BIN));

	// ---- Construct HitBuffer ----------------------------------------------
	SimpleThreadPool search_pool;
	Search::HitBuffer buf(
		key_partition, /*tmpdir=*/".", /*long_subject_offsets=*/false,
		QUERY_CONTEXTS, THREAD_COUNT,
		max_query, max_target, search_pool);

	// ---- Writing phase: one thread per bin --------------------------------
	std::atomic<uint64_t> expected_cs{ 0 };
	{
		std::vector<std::thread> writers;
		writers.reserve(THREAD_COUNT);
		for (int t = 0; t < THREAD_COUNT; ++t) {
			writers.emplace_back([&, t]() {
				Search::HitBuffer::Writer w(buf, (size_t)t);
				uint64_t local_cs = 0;
				const int q_begin = t * QUERIES_PER_BIN;
				const int q_end   = q_begin + QUERIES_PER_BIN;
				for (int q = q_begin; q < q_end; ++q) {
					const Loc      seed_off = (Loc)(q % 64);
					const uint32_t query    = (uint32_t)q;
					w.new_query(query, seed_off);
					for (int j = 0; j < HITS_PER_QUERY; ++j) {
						// subject in [1, TARGET_LEN]
						const uint32_t subj  = uint32_t((q * HITS_PER_QUERY + j) % TARGET_LEN + 1);
						// score in [1, 65534]
						const uint16_t score = uint16_t((q * 7 + j * 13) % 65534 + 1);
						w.write(query, PackedLoc(subj), score);
						local_cs += hit_fp(query, subj, (uint32_t)seed_off, score);
					}
				}
				expected_cs.fetch_add(local_cs, std::memory_order_relaxed);
				// Writer destructor flushes all bins automatically
			});
		}
		for (auto& t : writers)
			t.join();
	} // Writers destroyed here, flushing is complete

	buf.finish_writing();
	buf.alloc_buffer();

	// ---- Loading phase ----------------------------------------------------
	size_t   actual_total = 0;
	uint64_t actual_cs    = 0;

	while (buf.load(std::numeric_limits<size_t>::max())) {
		auto [hits, count, key_begin, key_end] = buf.retrieve();
		if (hits) {
			for (size_t i = 0; i < count; ++i) {
				const Search::Hit& h = hits[i];
				actual_cs += hit_fp(
					h.query_,
					(uint64_t)h.subject_,
					(uint32_t)h.seed_offset_,
					h.score_);
			}
			actual_total += count;
		}
	}

	buf.free_buffer();

	// ---- Restore config ---------------------------------------------------
	config.trace_pt_membuf = saved_membuf;
	config.swipe_all       = saved_swipe;

	const size_t   expected_total = (size_t)QUERY_COUNT * HITS_PER_QUERY;
	const uint64_t expected_sum   = expected_cs.load();
	const bool passed = (actual_total == expected_total) && (actual_cs == expected_sum);
	return { passed, expected_total, actual_total, expected_sum, actual_cs };
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

int run_hit_buffer_stress_test()
{
	std::cout << "\nHitBuffer stress test" << std::endl;
	std::cout << "=====================" << std::endl;
	std::cout << "Threads = " << config.threads_ << std::endl;

	int failures = 0;
	for (bool membuf : { true, false }) {
		const char* mode = membuf ? "in-memory (membuf)" : "disk";
		std::cout << "  Mode: " << mode << " ... " << std::flush;
		HitBufferTestResult r;
		try {
			r = run_single_mode(membuf);
		}
		catch (const std::exception& e) {
			std::cout << "EXCEPTION: " << e.what() << std::endl;
			++failures;
			continue;
		}
		if (r.passed) {
			std::cout << "PASSED"
				<< " (" << r.actual_hits << " hits, checksum ok)" << std::endl;
		}
		else {
			std::cout << "FAILED"
				<< "\n    expected hits=" << r.expected_hits
				<< " actual=" << r.actual_hits
				<< "\n    expected checksum=" << r.expected_checksum
				<< " actual=" << r.actual_checksum << std::endl;
			++failures;
		}
	}

	std::cout << "  Result: "
		<< (2 - failures) << "/2 passed" << std::endl;
	std::cout << "=====================" << std::endl;
	return failures;
}
