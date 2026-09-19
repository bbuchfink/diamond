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

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <stdint.h>
#include "lin_index.h"
#include "web_graph.h"
#include "basic/config.h"
#include "basic/seed.h"
#include "basic/shape_config.h"
#include "data/block/block.h"
#include "stats/score_matrix.h"
#include "search/hamming/finger_print.h"
#include "search/seed_array/enum_seeds.h"
#include "search/stage2.h"
#include "util/hash_function.h"
#include "util/log_stream.h"
#include "util/ptr_vector.h"
#include "util/simd/dispatch.h"

using std::endl;
using std::array;
using std::vector;
using ::DISPATCH_ARCH::FingerPrint;

namespace Search { namespace DISPATCH_ARCH {

/* Scans the member block against the seed index of the reference block. For
   every seed of the member block the index yields at most one location, the
   occurrence of that seed in the longest sequence of the reference block, which
   takes the role of the query in the resulting seed hit. */

/* Open addressing hash table with linear probing whose contents can be discarded in
   constant time: every slot carries the number of the scan it was written in, so that
   clearing the table is a single increment. Holds state that belongs to the member
   sequence currently being scanned and is thrown away when the next one begins. An
   entry carries a payload, which is an empty struct for the tables that only serve as
   a set. */
template<typename Payload>
struct EpochTable {

	EpochTable() :
		mask_((size_t(1) << INIT_BITS) - 1),
		count_(0),
		epoch_(1),
		table_(size_t(1) << INIT_BITS)
	{}

	void clear() {
		if (++epoch_ == 0) {
			// Wrapped around: no slot may keep a stamp that the new epoch could match.
			table_.assign(table_.size(), Slot());
			epoch_ = 1;
		}
		count_ = 0;
	}

	// Payload of the entry of the key, nullptr if the key is not contained. Invalidated
	// by the next insert.
	const Payload* find(const uint64_t key) const {
		for (size_t i = slot(key);; i = (i + 1) & mask_) {
			const Slot& s = table_[i];
			if (s.epoch != epoch_)
				return nullptr;
			if (s.key == key)
				return &s.payload;
		}
	}

	bool contains(const uint64_t key) const {
		return find(key) != nullptr;
	}

	// Inserts the key if it is not present and returns the payload of its entry.
	Payload& insert(const uint64_t key) {
		if (++count_ * 2 > (int64_t)table_.size())
			grow();
		for (size_t i = slot(key);; i = (i + 1) & mask_) {
			Slot& s = table_[i];
			if (s.epoch != epoch_) {
				s.epoch = epoch_;
				s.key = key;
				s.payload = Payload();
				return s.payload;
			}
			if (s.key == key)
				return s.payload;
		}
	}

private:

	struct Slot {
		uint64_t key = 0;
		uint32_t epoch = 0;
		Payload payload = Payload();
	};

	static constexpr int INIT_BITS = 10;

	size_t slot(const uint64_t key) const {
		return (size_t)(MurmurHash()(key) & (uint64_t)mask_);
	}

	void grow() {
		const std::vector<Slot> old(std::move(table_));
		table_.assign(old.size() * 2, Slot());
		mask_ = table_.size() - 1;
		for (const Slot& o : old) {
			if (o.epoch != epoch_)
				continue;
			size_t i = slot(o.key);
			while (table_[i].epoch == epoch_)
				i = (i + 1) & mask_;
			table_[i] = o;
		}
	}

	size_t mask_;
	int64_t count_;
	uint32_t epoch_;
	std::vector<Slot> table_;

};

struct EmptyPayload {};

/* Set of the diagonals on which a seed hit has been written for the member sequence
   that is currently being scanned. */
using DiagonalSet = EpochTable<EmptyPayload>;

/* Heuristic prefilter for the pair of sequences of a seed hit: the member sequence that
   is currently being scanned and the rep sequence that the index pointed at. It bounds
   the number of exactly matching 3-mers that an alignment passing the identity and
   coverage cutoffs of the round would have to contain, and rejects the pair if the two
   sequences do not share that many 3-mers at all.

   Let an accepted alignment have L columns of which I are identities.

   The coverage cutoffs put a lower bound l on L: the aligned range of either sequence
   spans at most L columns, so --query-cover bounds L below by its fraction of the rep
   length, --subject-cover by its fraction of the member length, and the cutoff of
   --member-cover, which the round passes to the search as --query-or-target-cover, by
   its fraction of the shorter of the two.

   The identity cutoff gives I >= id * L. Every column that is not an identity ends a run
   of identities, so the I identities form at most L - I + 1 runs, and a run of r
   identities contains r - 2 exactly matching 3-mers. The alignment therefore contains at
   least

       I - 2 * (L - I + 1) = (3 * id - 2) * L - 2 >= (3 * id - 2) * l - 2

   of them, the coefficient being positive because the filter is only used from MIN_ID
   on. Distinct matching 3-mers of an alignment occupy distinct positions in either
   sequence, so their number cannot exceed the size of the multiset intersection of the
   3-mers of the two sequences, which is what is counted here.

   The count ignores where in the sequences the 3-mers sit and in which order, which is
   what makes it cheap, and the bound is a heuristic rather than an exact one: the
   identity the round cuts on is estimated from the score (Stats::approx_id) and is not
   the fraction of identical columns, and gaps inside the aligned ranges let L exceed the
   coverage bound in ways the argument above does not see. Some true pairs are lost.

   The 3-mers of the member sequence are counted in a table addressed by the 3-mer
   itself, five bits per letter, and the rep sequence is scanned against it. The member
   sequence is the same for a whole run of seed hits, so the table is built once per
   member sequence, and the verdict for a rep sequence is cached for that run. Low
   complexity regions are hard masked in the input of a clustering run; their letters are
   counted like any others, which can only inflate the number of shared 3-mers and weaken
   the filter, never reject a pair that would otherwise be accepted.

   How aggressive the filter is can be raised with --kmer-prefilter a, at the cost of more
   lost pairs. The coefficient of the bound is moved from the worst case 3 * id - 2, in
   which the mismatches are spread as evenly as possible, towards id^3, the expected
   density of matching 3-mers if the mismatches fall at random:

       coefficient = (3 * id - 2) + a * (id^3 - (3 * id - 2))

   a = 0 is the bound above and the default. At a = 1 the filter asks for as many 3-mers
   as a typical alignment exactly at the cutoffs contains, so that a good part of the
   pairs close to the cutoffs is lost; values above 1 go beyond that. Interpolating
   rather than scaling the bound keeps the meaning of a comparable across identity
   cutoffs: the gap between the two coefficients is 0.24 at an identity of 70% but only
   0.03 at 90%. A negative value turns the filter off. */
struct KmerFilter {

	enum class Verdict : uint8_t { UNKNOWN, PASS, FAIL };

	/* Identity cutoff below which the filter is not used: the coefficient 3 * id - 2 of
	   the bound vanishes at 2/3, and nowhere near it is the bound expected to be of any
	   use. */
	static constexpr double MIN_ID = 70.0;
	static constexpr int K = 3;
	static constexpr int LETTER_BITS = 5;
	static constexpr uint32_t KMER_MASK = (uint32_t(1) << (K * LETTER_BITS)) - 1;

	/* Not used with --mutual-cover, with a negative --kmer-prefilter, or below MIN_ID.
	   Without any coverage cutoff there is no lower bound on the length of the alignment
	   and nothing to filter on. */
	static bool enabled() {
		if (config.kmer_prefilter < 0.0)
			return false;
		if (config.mutual_cover.present())
			return false;
		if (id_cutoff() < MIN_ID)
			return false;
		return std::max({ config.query_cover, config.subject_cover, config.query_or_target_cover }) > 0.0;
	}

	KmerFilter() :
		enabled_(enabled()),
		id_coeff_(coefficient(id_cutoff() / 100.0, config.kmer_prefilter)),
		query_cover_(config.query_cover / 100.0),
		subject_cover_(config.subject_cover / 100.0),
		either_cover_(config.query_or_target_cover / 100.0),
		member_len_(0),
		test_count_(0)
	{
		if (enabled_)
			count_.resize((size_t)KMER_MASK + 1, 0);
	}

	// Begins the run of seed hits of a member sequence: counts its 3-mers and discards
	// the verdicts of the previous member sequence.
	void set_member(const Sequence& seq) {
		if (!enabled_)
			return;
		for (const uint32_t k : member_kmers_)
			count_[k] = 0;
		member_kmers_.clear();
		cache_.clear();
		member_len_ = seq.length();
		uint32_t kmer = 0;
		for (Loc i = 0; i < member_len_; ++i) {
			kmer = ((kmer << LETTER_BITS) | (uint32_t)seq[i]) & KMER_MASK;
			if (i < K - 1)
				continue;
			uint16_t& c = count_[kmer];
			if (c == 0)
				member_kmers_.push_back(kmer);
			if (c < COUNT_MAX)
				++c;
		}
	}

	/* Verdict for the rep sequence if it has been computed for this member sequence
	   before. Lets the remaining seed hits of a pair that has already been rejected be
	   discarded without running the ungapped filters on them. */
	Verdict cached(const unsigned rep_id) const {
		if (!enabled_)
			return Verdict::UNKNOWN;
		const Verdict* v = cache_.find(rep_id);
		return v == nullptr ? Verdict::UNKNOWN : *v;
	}

	bool test(const unsigned rep_id, const Sequence& rep) {
		if (!enabled_)
			return true;
		if (const Verdict* v = cache_.find(rep_id))
			return *v == Verdict::PASS;
		++test_count_;
		const int64_t required = required_kmers(rep.length());
		const bool pass = required <= 0 || shared_kmers(rep) >= required;
		cache_.insert(rep_id) = pass ? Verdict::PASS : Verdict::FAIL;
		return pass;
	}

	int64_t test_count() const {
		return test_count_;
	}

private:

	// Saturation value of the 3-mer counts of the member sequence. Saturating can only
	// lower the bound of a pair, never reject one that would otherwise be accepted.
	static constexpr uint16_t COUNT_MAX = 0xffff;

	static double id_cutoff() {
		return std::max(config.approx_min_id.get(0.0), config.min_id);
	}

	// Coefficient of the bound for an identity (as a fraction) and the aggressiveness a
	// of --kmer-prefilter, see above.
	static double coefficient(const double id, const double a) {
		const double worst_case = 3.0 * id - 2.0;
		return worst_case + std::max(a, 0.0) * (id * id * id - worst_case);
	}

	// Number of exactly matching 3-mers that an accepted alignment against a rep
	// sequence of that length has to contain, see above.
	int64_t required_kmers(const Loc rep_len) const {
		const double l = std::max({ query_cover_ * rep_len, subject_cover_ * member_len_,
			either_cover_ * std::min(rep_len, member_len_) });
		return (int64_t)(id_coeff_ * l) - (K - 1);
	}

	/* Size of the multiset intersection of the 3-mers of the rep sequence and those of
	   the member sequence, i.e. the number of 3-mers of the rep sequence that can be
	   paired with a distinct occurrence of the same 3-mer in the member sequence. The
	   counts of the member table are consumed while the rep sequence is walked and
	   restored afterwards. */
	int64_t shared_kmers(const Sequence& rep) {
		int64_t n = 0;
		uint32_t kmer = 0;
		const Loc len = rep.length();
		for (Loc i = 0; i < len; ++i) {
			kmer = ((kmer << LETTER_BITS) | (uint32_t)rep[i]) & KMER_MASK;
			if (i < K - 1)
				continue;
			uint16_t& c = count_[kmer];
			if (c > 0) {
				--c;
				++n;
				taken_.push_back(kmer);
			}
		}
		for (const uint32_t k : taken_)
			++count_[k];
		taken_.clear();
		return n;
	}

	const bool enabled_;
	const double id_coeff_, query_cover_, subject_cover_, either_cover_;
	Loc member_len_;
	int64_t test_count_;
	std::vector<uint16_t> count_;
	std::vector<uint32_t> member_kmers_, taken_;
	EpochTable<Verdict> cache_;

};

/* Vertices of the graph of the pairs of sequences that have been paired by a seed hit,
   see scan_lin_index. The rep sequences are numbered first, followed by the member
   sequences, unless both roles are taken by the same block, in which case a sequence is
   the same vertex in either role and a pair is recognized in both orientations. */
struct PairVertices {

	PairVertices(const Search::Config& cfg) :
		self(config.self && cfg.current_ref_block == 0),
		rep_count((WebGraph::Vertex)cfg.query->seqs().size()),
		member_count((WebGraph::Vertex)cfg.target->seqs().size())
	{}

	WebGraph::Vertex count() const {
		return self ? std::max(rep_count, member_count) : rep_count + member_count;
	}

	WebGraph::Vertex rep(const unsigned rep_id) const {
		return (WebGraph::Vertex)rep_id;
	}

	WebGraph::Vertex member(const uint32_t member_id) const {
		return self ? (WebGraph::Vertex)member_id : rep_count + (WebGraph::Vertex)member_id;
	}

	const bool self;
	const WebGraph::Vertex rep_count, member_count;

};

struct ScanCallback {

	/* known_pairs holds the pairs of sequences that the previous seed shapes have written
	   seed hits for, nullptr if there are none. With record_pairs set, the pairs that
	   this shape writes seed hits for are collected in new_pairs. */
	ScanCallback(const Context& context, Search::Config& cfg, const LinIndex& index, unsigned shape_id, size_t thread_id,
		const WebGraph* known_pairs, bool record_pairs) :
		writer(*cfg.seed_hit_buf, thread_id),
		work_set(context, cfg, shape_id, &writer, nullptr, nullptr),
		index(index),
		query_seqs(cfg.query->seqs()),
		target_seqs(cfg.target->seqs()),
		target_seed_hits(cfg.target_seed_hits.get()),
		hamming_filter_id(cfg.hamming_filter_id),
		self(config.self && cfg.current_ref_block == 0),
		vertices(cfg),
		known_pairs(known_pairs),
		record_pairs(record_pairs),
		lookup_count(0),
		paired_count(0),
		dup_diag_count(0),
		kmer_filtered_count(0),
		known_pair_count(0)
	{}

	/* Diagonals of different sequences of the reference block are unrelated, so the
	   sequence is part of the key. The member sequence is the same for all keys of
	   one scan and does not have to be included. */
	static uint64_t diagonal_key(const unsigned query_id, const Loc diagonal) {
		return ((uint64_t)query_id << 32) | (uint64_t)(uint32_t)diagonal;
	}

	bool operator()(const uint64_t key, const uint64_t pos, const uint32_t block_id, const uint64_t shape_id) {
		++lookup_count;
		const uint64_t pivot = index(key);
		if (pivot == 0)
			return true;
		work_set.stats.inc(Statistics::SEED_HITS);

		const std::pair<BlockId, Loc> l = query_seqs.local_position((int64_t)pivot);
		const unsigned query_id = l.first;
		if (self && query_id == block_id)
			return true;

		++paired_count;

		/* The seeds are enumerated one member sequence after the other, so a change of
		   block_id begins a new set of diagonals. */
		if (block_id != last_block_id_) {
			diagonals.clear();
			recorded_reps.clear();
			member_vertex_ = vertices.member(block_id);
			kmer_filter.set_member(target_seqs[block_id]);
			member_seq_begin_ = target_seqs.position(block_id, 0);
			last_block_id_ = block_id;
		}

		/* The pair of sequences has been written by a previous seed shape, on whichever
		   diagonal: it is known to the search already and is not looked at again. The
		   seed counts as hit, like one on a redundant diagonal below. */
		if (known_pairs && known_pairs->contains(vertices.rep(query_id), member_vertex_)) {
			++known_pair_count;
			if (target_seed_hits)
				target_seed_hits->operator[]((size_t)shape_id).atomic_set(pos);
			return true;
		}

		const Loc seed_offset = l.second;
		/* A seed hit on a diagonal that has already produced a hit is redundant: it
		   extends into the same ungapped alignment. The diagonal is recorded when the
		   hit is written, not here, so that a diagonal whose first seed is rejected by
		   the filters below can still be found by a later seed. */
		const uint64_t diagonal = diagonal_key(query_id, seed_offset - Loc((int64_t)pos - member_seq_begin_));
		if (diagonals.contains(diagonal)) {
			++dup_diag_count;
			if (target_seed_hits)
				target_seed_hits->operator[]((size_t)shape_id).atomic_set(pos);
			return true;
		}

		// The pair of sequences has already been rejected by the 3-mer filter.
		if (kmer_filter.cached(query_id) == KmerFilter::Verdict::FAIL) {
			++kmer_filtered_count;
			return true;
		}

		alignas(64) array<char, 48> fq, fs;
		FingerPrint::load(query_seqs.data(pivot), &fq);
		FingerPrint::load(target_seqs.data(pos), &fs);
		if (FingerPrint(fq).match(FingerPrint(fs)) < hamming_filter_id)
			return true;
		work_set.stats.inc(Statistics::TENTATIVE_MATCHES1);

		const int query_len = query_seqs.length(query_id);
		const int score_cutoff = ungapped_cutoff(query_len, work_set);
		int score = std::numeric_limits<int>::max();
		if (score_cutoff) {
			const int window = ungapped_window(query_len);
			const Letter* query = query_seqs.data(pivot);
			const Sequence query_clipped = Util::Seq::clip(query - window, window * 2, window);
			const int window_left = int(query - query_clipped.data()), window_clipped = (int)query_clipped.length();
			const Letter* subject = target_seqs.data(pos) - window_left;
			DP::window_ungapped_best(query_clipped.data(), &subject, 1, window_clipped, &score);
			if (score <= score_cutoff)
				return true;
		}
		work_set.stats.inc(Statistics::TENTATIVE_MATCHES3);

		/* The 3-mer filter is linear in the length of the rep sequence and therefore the
		   most expensive of the filters, so it runs last. Its verdict holds for the whole
		   pair of sequences and is computed at most once per rep sequence and member
		   sequence. */
		if (!kmer_filter.test(query_id, query_seqs[query_id])) {
			++kmer_filtered_count;
			return true;
		}

		work_set.stats.inc(Statistics::TENTATIVE_MATCHES3B);

		if (target_seed_hits)
			target_seed_hits->operator[]((size_t)shape_id).atomic_set(pos);
		diagonals.insert(diagonal);
		if (record_pairs && !recorded_reps.contains(query_id)) {
			recorded_reps.insert(query_id);
			new_pairs.emplace_back(vertices.rep(query_id), member_vertex_);
		}
		if (query_id != last_query_ || seed_offset != last_seed_offset_) {
			writer.new_query(query_id, seed_offset);
			last_query_ = query_id;
			last_seed_offset_ = seed_offset;
		}
		writer.write(query_id, PackedLoc(pos), (uint16_t)score, block_id);
		return true;
	}

	void finish() {}

	HitBuffer::Writer writer;
	WorkSet work_set;
	const LinIndex& index;
	const SequenceSet& query_seqs, & target_seqs;
	std::vector<BitVector>* const target_seed_hits;
	const unsigned hamming_filter_id;
	const bool self;
	const PairVertices vertices;
	const WebGraph* const known_pairs;
	const bool record_pairs;
	DiagonalSet diagonals;
	// Rep sequences whose pair with the current member sequence is in new_pairs.
	DiagonalSet recorded_reps;
	std::vector<WebGraph::Edge> new_pairs;
	KmerFilter kmer_filter;
	unsigned last_query_ = std::numeric_limits<unsigned>::max();
	Loc last_seed_offset_ = std::numeric_limits<Loc>::min();
	uint32_t last_block_id_ = std::numeric_limits<uint32_t>::max();
	int64_t member_seq_begin_ = 0;
	WebGraph::Vertex member_vertex_ = 0;
	uint64_t lookup_count, paired_count, dup_diag_count, kmer_filtered_count, known_pair_count;

};

/* The seed shapes are processed one after the other: the hash table of the reference
   block is rebuilt for the shape, the member block is scanned for that shape only, and
   the table is released again before the next shape is loaded. Holding all shapes at
   once would multiply the memory of the index by their number.

   The pairs of a rep and a member sequence that a shape writes seed hits for are
   collected in an undirected graph that is carried over to the following shapes. A seed
   hit of a later shape for a pair that is in the graph is discarded before any of the
   filters run, irrespective of its diagonal: the pair has been handed to the search
   already. Within a shape, a pair can still produce seed hits on several diagonals. The
   graph is only read while a shape is scanned and grows in between, so that the threads
   need no synchronization. */
void scan_lin_index(Search::Config& cfg) {
	const int threads = std::max(config.threads_, 1);
	cfg.lin_index->check_encoding(cfg);
	const auto partition = cfg.target->seqs().partition(threads);
	WebGraph known_pairs(PairVertices(cfg).count());

	for (int s = 0; s < shapes.count(); ++s) {
		TaskTimer timer;
		const std::string label = " (shape " + std::to_string(s) + ")";

		/* Soft masking is applied to the reference block when its index is built. The
		   member block is scanned unmasked so that the fingerprints of both blocks are
		   computed on the same, unmasked sequence data. The seed encoding is the one the
		   index was built with, see LinIndex::encoding; filter_masked_seeds makes the
		   hashed encoding reject a window with a masked letter on a match position of the
		   shape, which is what the letterwise encoding does unconditionally. */
		const EnumCfg enum_cfg{ &partition, s, s + 1, cfg.lin_index->seed_encoding(), nullptr, true, false, cfg.seed_complexity_cut,
			MaskingAlgo::NONE, cfg.minimizer_window, false, false, cfg.sketch_size, cfg.target_seed_hits.get() };

		cfg.lin_index->build_shape(s);

		*message_stream << "Rep block: " << cfg.query->seqs().mem_size() << ", member block: " << cfg.target->seqs().mem_size() << ", index size: "
			<< cfg.lin_index->table_size() << ", slots used: " << (int64_t)cfg.lin_index->insert_count() * LinIndex::Table::SLOT_BYTES << endl;
		const vector<uint32_t> patterns = shapes.patterns(0, s + 1);
		// A Context holds two PatternMatcher lookup tables of 2^MAX_SHAPE_LEN bytes
		// each and must not be placed on the stack.
		const std::unique_ptr<Context> context(new Context{ { patterns.data(), patterns.data() + patterns.size() - 1 },
			{ patterns.data(), patterns.data() + patterns.size() },
			score_matrix.rawscore(config.short_query_ungapped_bitscore),
			nullptr,
			seedp_mask(cfg.seedp_bits) });

		timer.go("Scanning member block" + label);
		PtrVector<ScanCallback> v;
		for (int i = 0; i < threads; ++i)
			v.push_back(new ScanCallback(*context, cfg, *cfg.lin_index, (unsigned)s, i, s > 0 ? &known_pairs : nullptr, s + 1 < shapes.count()));
		enum_seeds(*cfg.target, v, &no_filter, enum_cfg);
		uint64_t lookup_count = 0, paired_count = 0, dup_diag_count = 0, kmer_filtered_count = 0, known_pair_count = 0;
		int64_t kmer_test_count = 0, pair_buffer_bytes = 0;
		vector<WebGraph::Edge> new_pairs;
		for (int i = 0; i < threads; ++i) {
			statistics += v[i].work_set.stats;
			lookup_count += v[i].lookup_count;
			paired_count += v[i].paired_count;
			dup_diag_count += v[i].dup_diag_count;
			kmer_filtered_count += v[i].kmer_filtered_count;
			known_pair_count += v[i].known_pair_count;
			kmer_test_count += v[i].kmer_filter.test_count();
			pair_buffer_bytes += (int64_t)(v[i].new_pairs.size() * sizeof(WebGraph::Edge));
			new_pairs.insert(new_pairs.end(), v[i].new_pairs.begin(), v[i].new_pairs.end());
		}
		v.clear();
		timer.finish();
		*message_stream << "Hit seeds: " << paired_count << '/' << lookup_count << " (" << (lookup_count ? (100.0 * paired_count) / lookup_count : 0)
			<< "%), redundant diagonals: " << dup_diag_count << ", known pairs: " << known_pair_count << std::endl;
		if (KmerFilter::enabled())
			*message_stream << "Seed hits filtered by the 3-mer bound: " << kmer_filtered_count << ", sequence pairs tested: " << kmer_test_count << endl;
		if (!new_pairs.empty()) {
			timer.go("Merging sequence pairs" + label);
			known_pairs.insert(new_pairs);
			timer.finish();
		}
		*log_stream << "Pair graph" << label << ": edges=" << known_pairs.edge_count() << " size=" << known_pairs.data_size()
			<< " bytes, per-thread pair buffers=" << pair_buffer_bytes << " bytes" << endl;
	}
	cfg.lin_index->free_shape();
}

}

DISPATCH_1V(scan_lin_index, Search::Config&, cfg)

}
