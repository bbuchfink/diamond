/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#include <inttypes.h>
#include <algorithm>
#include <memory>
#include <numeric>
#include <fstream>
#include <sstream>
#include "multinode.h"
#include "search/lin_index/lin_index.h"
#include "search/seed_array/enum_seeds.h"
#include "multinode.h"
#include "util/algo/hyperloglog.h"

using std::endl;
using std::string;
using std::ofstream;
using std::string;
using std::unique_ptr;
using std::unordered_map;
using std::runtime_error;
using std::vector;
using std::pair;

// Name of the file holding the seed counter of the superblocks, written next to their
// volume list.
static const char* SEED_COUNT_FILE = "seed_counts.tsv";

struct Cost {
	uint64_t mem_use, io;
};

/* Memory taken by the seed index of one superblock, given the number of distinct seeds
   of its most frequent shape (only one shape is resident at a time). The index is far
   larger while it is built than while it is resident for the scan of a member block: the
   build table keeps the sequence lengths in an extra two bytes per slot and is sized from
   an estimate plus a safety margin, so the two differ by a factor of 1.5625. Both are
   needed because they occur at different points of the schedule. */
struct IndexSize {
	uint64_t build, search;
};

static IndexSize index_size(const uint64_t distinct_seeds) {
	return { Search::LinIndex::build_memory(distinct_seeds), Search::LinIndex::index_memory(distinct_seeds) };
}

/* Memory taken by a superblock in its role as the incoming member block, i.e. the one
   that is scanned against the seed indices of all the superblocks before it and does not
   need an index of its own. Besides its letters it holds the seed hit bit vectors
   (`target_seed_hits`), one bit per letter and seed shape, which live for as long as the
   block is the member. */
static uint64_t member_block_memory(const uint64_t letters) {
	return letters + (letters + 63) / 64 * 8 * shapes.count();
}

/* Estimated number of distinct seeds of one superblock, maximized over the seed shapes,
   obtained by merging the counters of the minichunks the superblock is made of. It is the
   number of pivots that the seed index of the superblock will hold for its most frequent
   shape, and all of its hash tables are sized from it: only one shape is resident at a
   time, so nothing is gained by sizing them one by one. */
static uint64_t superblock_seed_count(vector<vector<HyperLogLog>>::const_iterator begin, vector<vector<HyperLogLog>>::const_iterator end) {
	if (begin == end)
		return 0;
	uint64_t count = 0;
	for (int i = 0; i < shapes.count(); ++i) {
		HyperLogLog m(begin->front().precision());
		for (auto it = begin; it != end; ++it)
			m.merge(it->operator[](i));
		count = std::max(count, (uint64_t)m.estimate());
	}
	return count;
}

static IndexSize mem_use_index(vector<vector<HyperLogLog>>::const_iterator begin, vector<vector<HyperLogLog>>::const_iterator end) {
	return index_size(superblock_seed_count(begin, end));
}

string seed_count_file(const string& volume_list_file) {
	return containing_directory(volume_list_file) + PATH_SEPARATOR + SEED_COUNT_FILE;
}

unordered_map<string, uint64_t> read_seed_counts(const string& file_name) {
	unordered_map<string, uint64_t> counts;
	std::ifstream in(file_name);
	if (!in.good())
		return counts;
	string line, path;
	uint64_t n;
	while (std::getline(in, line)) {
		std::istringstream row(line);
		if (row >> path >> n)
			counts[path] = n;
	}
	return counts;
}

// letters_a belongs to the superblock whose index is resident, letters_b to the member block.
static uint64_t mem_use_combo(uint64_t index_size, uint64_t letters_a, uint64_t letters_b) {
	return index_size + letters_a + member_block_memory(letters_b);
}

static Cost cost_member_superblock(size_t member_superblock, const vector<uint64_t>& letters, const vector<IndexSize>& index_sizes) {
	uint64_t total = 0, letters_member = letters[member_superblock], io = letters_member;
	for (size_t i = 0; i < member_superblock; ++i) {
		total = std::max(total, mem_use_combo(index_sizes[i].search, letters[i], letters_member));
		io += letters[i];
	}
	// The superblock is loaded on its own while its index is built, which is the point at
	// which the larger of the two index sizes is resident.
	return { std::max(total, letters_member + index_sizes[member_superblock].build), io };
}

static Cost cost_partition(const vector<size_t>& partition, const vector<vector<HyperLogLog>>& counters, const vector<uint64_t>& letters) {
	vector<IndexSize> index_sizes;
	vector<uint64_t> superblock_letters;
	for (size_t i = 0; i < partition.size() - 1; ++i) {
		index_sizes.push_back(mem_use_index(counters.begin() + partition[i], counters.begin() + partition[i + 1]));
		superblock_letters.push_back(0);
		for (size_t j = partition[i]; j < partition[i + 1]; ++j)
			superblock_letters.back() += letters[j];
	}
	uint64_t total = 0, io = 0;
	for (size_t i = 0; i < partition.size() - 1; ++i) {
		const Cost c = cost_member_superblock(i, superblock_letters, index_sizes);
		total = std::max(total, c.mem_use);
		io += c.io;
	}
	return { total, io };
}

/* Partitioning the minichunks into superblocks.

The cost model above decomposes exactly. Write L_m for the letters of superblock
m, B_m = L_m + build index size(m) for the peak while the index of m is built,
M_m = member_block_memory(L_m) for what m takes while it is the member block,
S_m = L_m + search index size(m) for what stays resident while m is the rep of a
later member, prefS_m = max over j < m of S_j (0 for m = 0), and let cum[] be the
prefix sums of the per minichunk letter counts. Then

	io(P)  = cum[n] + sum of cum[q] over the internal cut positions q of P
	mem(P) = max over m of max(B_m, prefS_m + M_m)

So the io cost is just the sum of the cumulative letter counts at the cut points:
every cut is paid for independently, and late cuts are far more expensive than
early ones. The memory cost couples the superblocks only through the running
maximum of S, while B is local to a superblock. Both facts together allow a left
to right dynamic program over the cut positions whose state is (running max of S,
io so far). Refining a partition can only lower mem() and raise io(), so the
search is a walk along that trade off.

The state set at each position is reduced to its Pareto front (both components
are minimized) and then quantized to at most BEAM entries, which is what makes
this a heuristic rather than an exact solver. The state with the smallest running
max of A is always kept, so a feasible partition is found whenever one exists.
*/

static const size_t BEAM = 64;
// A state set is pruned on arrival, but it accumulates states from every earlier position
// that can reach it, so it is also pruned whenever it grows past this. That bounds the
// memory of the dp by O(n * PRUNE_AT) states instead of O(n * BEAM * superblock size).
static const size_t PRUNE_AT = 8 * BEAM;

struct DpState {
	uint64_t pref_a;                // max S over the superblocks chosen so far
	uint64_t io;                    // sum of cum[] over the cuts made so far
	uint64_t mem;                   // memory cost of the partial partition, for tie breaking only
	size_t prev_pos, prev_state;    // backtracking
};

static void pareto_front(vector<DpState>& v, const uint64_t memory_limit) {
	if (v.size() <= 1)
		return;
	std::sort(v.begin(), v.end(), [](const DpState& a, const DpState& b) {
		return a.pref_a != b.pref_a ? a.pref_a < b.pref_a : (a.io != b.io ? a.io < b.io : a.mem < b.mem); });
	vector<DpState> f;
	for (const DpState& s : v)
		if (f.empty() || s.io < f.back().io)
			f.push_back(s);
	if (f.size() > BEAM) {
		// The only thing about a state that matters for the future is its running max of S,
		// through the capacity memory_limit - pref_a it leaves for the superblocks still to
		// come. So quantize pref_a into BEAM bands and keep the smallest io per band: a kept
		// state then stands in for a discarded one whose pref_a is at most one band width
		// smaller. f[0] is pinned because it has the smallest pref_a of all, which is what
		// guarantees that a feasible partition is still found whenever one exists.
		const uint64_t width = memory_limit / BEAM + 1;
		vector<DpState> g;
		g.push_back(f[0]);
		for (size_t i = 1; i < f.size(); ++i)
			if (g.size() > 1 && g.back().pref_a / width == f[i].pref_a / width)
				g.back() = f[i];
			else
				g.push_back(f[i]);
		f = std::move(g);
	}
	v = std::move(f);
}

/* B(begin, r) and S(begin, r) = index size + letters of the superblock [begin, r), for
   r = begin + 1 onwards, truncated as soon as B exceeds the memory limit (B >= S always,
   so B is what decides whether the superblock is feasible on its own). Computed
   incrementally, so building all rows costs one HyperLogLog merge per feasible superblock
   instead of one per member. */
struct SuperblockCost {
	uint64_t build, search;
};

static vector<SuperblockCost> superblock_cost_row(const size_t begin, const vector<vector<HyperLogLog>>& counters, const vector<uint64_t>& cum, const uint64_t memory_limit) {
	const size_t n = counters.size();
	vector<SuperblockCost> row;
	vector<HyperLogLog> acc = counters[begin];
	for (size_t r = begin + 1; r <= n; ++r) {
		if (r > begin + 1)
			for (size_t i = 0; i < acc.size(); ++i)
				acc[i].merge(counters[r - 1][i]);
		uint64_t distinct = 0;
		for (const HyperLogLog& h : acc)
			distinct = std::max(distinct, (uint64_t)h.estimate());
		const IndexSize idx = index_size(distinct);
		const uint64_t letters = cum[r] - cum[begin];
		if (idx.build + letters > memory_limit)
			break;
		row.push_back({ idx.build + letters, idx.search + letters });
	}
	return row;
}

static vector<size_t> combine_blocks(Job& job, const uint64_t memory_limit, const vector<vector<HyperLogLog>>& counters, const vector<uint64_t>& letters) {
	const size_t n = counters.size();
	if (n <= 1)
		return { 0, n };

	vector<uint64_t> cum(n + 1, 0);
	for (size_t i = 0; i < n; ++i)
		cum[i + 1] = cum[i] + letters[i];

	vector<vector<DpState>> front(n + 1);
	front[0].push_back({ 0, 0, 0, 0, 0 });
	for (size_t q = 0; q < n; ++q) {
		if (front[q].empty())
			continue;
		pareto_front(front[q], memory_limit);
		const vector<SuperblockCost> a = superblock_cost_row(q, counters, cum, memory_limit);
		for (size_t i = 0; i < front[q].size(); ++i) {
			const DpState s = front[q][i];
			for (size_t j = 0; j < a.size(); ++j) {
				const size_t r = q + 1 + j;
				const uint64_t l = cum[r] - cum[q];
				// The first superblock is never a member, so it pays for neither a resident
				// rep index nor the seed hit bit vectors.
				const uint64_t member = q > 0 ? s.pref_a + member_block_memory(l) : 0;
				if (member > memory_limit)
					break; // member is increasing in r
				front[r].push_back({ std::max(s.pref_a, a[j].search),
					s.io + (r < n ? cum[r] : 0),
					std::max(s.mem, std::max(a[j].build, member)),
					q, i });
				if (front[r].size() >= PRUNE_AT)
					pareto_front(front[r], memory_limit); // r > q, never the set being expanded
			}
		}
	}

	if (front[n].empty()) {
		// Not even one superblock per minichunk fits, and that is the memory minimal partition.
		job.log("Warning: no partition of the minichunks satisfies the memory limit of %" PRIu64 " bytes.", memory_limit);
		vector<size_t> p(n + 1);
		std::iota(p.begin(), p.end(), (size_t)0);
		return p;
	}

	size_t best = 0;
	for (size_t i = 1; i < front[n].size(); ++i)
		if (front[n][i].io < front[n][best].io || (front[n][i].io == front[n][best].io && front[n][i].mem < front[n][best].mem))
			best = i;
	vector<size_t> p;
	for (size_t pos = n, i = best;;) {
		p.push_back(pos);
		if (pos == 0)
			break;
		const DpState& s = front[pos][i];
		pos = s.prev_pos;
		i = s.prev_state;
	}
	std::reverse(p.begin(), p.end());

	const Cost c = cost_partition(p, counters, letters);
	job.log("Combining %" PRIu64 " minichunks into %" PRIu64 " superblocks. Memory=%" PRIu64 " limit=%" PRIu64 " IO=%" PRIu64,
		(uint64_t)n, (uint64_t)(p.size() - 1), c.mem_use, memory_limit, c.io);
	job.log("Minichunk partition: ");
	for (size_t i = 0; i < p.size() - 1; ++i)
		job.log("  %lli", (long long)p[i + 1] - (long long)p[i]);
	return p;
}

string make_merged_blocks(Job& job, const string& minichunks, const string& merged_dir, uint64_t letter_count) {
	static const double LETTER_OVERCOUNT = 0.256;
	const string vols_file = merged_dir + "volumes.tsv";
	Atomic combine_lock(job.base_dir() + "combine_lock", job), combine_done(job.base_dir() + "combine_done", job);
	if (combine_lock.fetch_add() == 0 && combine_done.get() == 0) {
		VolumedFile vols_in(minichunks);
		configure_round(job, 0);
		mkdir(merged_dir);
		ofstream vols(vols_file);
		if (!job.is_linear_round() || config.mutual_cover.present()) {
			if (!config.mutual_cover.present() && vols_in.size() > 1)
				throw runtime_error(">1 minichunk in non-linear round.");
			for (const Volume& v : vols_in)
				vols << v.path << endl;
			combine_done.fetch_add();
			return vols_file;
		}
		vector<vector<HyperLogLog>> counters;
		vector<uint64_t> letters;
		std::tie(counters, letters) = count_distinct_seeds(job, vols_in, letter_count);
		for (uint64_t& l : letters)
			l += l * LETTER_OVERCOUNT;
		const vector<size_t> superblocks = combine_blocks(job, job.mem_limit, counters, letters);
		
		// One line per superblock: its path followed by the estimated number of distinct
		// seeds of its most frequent shape, which is what the hash tables of its seed index
		// have to hold.
		ofstream counts(merged_dir + SEED_COUNT_FILE);
		for (size_t p = 0; p < superblocks.size() - 1; ++p) {
			const string path = merged_dir + "superblock" + std::to_string(p) + ".tsv";
			ofstream vol1(path);
			uint64_t n = 0;
			for (uint64_t i = superblocks[p]; i < superblocks[p + 1]; ++i) {
				vol1 << vols_in[i].path << endl; // << '\t' << 0 << endl;
				n += 0;
			}
			vols << path << '\t' << n << endl;
			counts << path << '\t' << superblock_seed_count(counters.begin() + superblocks[p], counters.begin() + superblocks[p + 1]) << endl;
		}

		combine_done.fetch_add();
	}
	else
		combine_done.await(1);	
	return vols_file;
}