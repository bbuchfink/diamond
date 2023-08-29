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

#include <array>
#include <thread>
#include <atomic>
#include "seed_complexity.h"
#include "../data/block/block.h"
#include "../util/algo/join_result.h"
#include "../util/string/string.h"

using std::array;
using std::vector;
using std::endl;

extern double lnfact[];
std::unordered_set<uint64_t> soft_mask;

bool Search::seed_is_complex(const Letter* seq, const Shape& shape, const double cut)
{
	array<unsigned, TRUE_AA> count;
	count.fill(0);
	for (int i = 0; i < shape.weight_; ++i) {
		const Letter l = letter_mask(seq[shape.positions_[i]]);
		if (l >= TRUE_AA)
			return false;
		++count[Reduction::reduction(l)];
	}
	double entropy = lnfact[shape.weight_];
	for (unsigned i = 0; i < Reduction::reduction.size(); ++i)
		entropy -= lnfact[count[i]];
	return entropy >= cut;
}

bool Search::seed_is_complex_unreduced(Letter* seq, const Shape& shape, const double cut, const bool mask_seeds, SeedStats& stats)
{
	array<unsigned, TRUE_AA> count;
	count.fill(0);
	for (int i = 0; i < shape.weight_; ++i) {
		const Letter l = letter_mask(seq[shape.positions_[i]]);
		if (l >= TRUE_AA) {
			if (mask_seeds) *seq |= SEED_MASK;
			return false;
		}
		++count[(int)l];
	}
	++stats.good_seed_positions;
	double entropy = lnfact[shape.weight_];
	for (unsigned i = 0; i < TRUE_AA; ++i)
		entropy -= lnfact[count[i]];
	if (entropy < cut) {
		if (mask_seeds) *seq |= SEED_MASK;
		++stats.low_complexity_seeds;
		return false;
	}
	return true;
}

void Search::mask_seeds(const Shape& shape, const SeedPartitionRange& range, DoubleArray<SeedLoc>* query_seed_hits, DoubleArray<SeedLoc>* ref_seed_hits, Search::Config& cfg)
{
	if (cfg.seed_encoding != SeedEncoding::SPACED_FACTOR)
		return;

	TaskTimer timer("Masking low complexity seeds");
	SequenceSet& query_seqs = cfg.query->seqs();
	std::atomic_int seedp(range.begin());
	std::atomic_size_t seed_count(0), masked_seed_count(0), query_count(0), target_count(0);
	const double cut = cfg.seed_complexity_cut;

	auto worker = [&] {
		int p;
		size_t sc(0), msc(0), qc(0), tc(0);
		while ((p = seedp++) < range.end()) {
			for (auto it = JoinIterator<SeedLoc>(query_seed_hits[p].begin(), ref_seed_hits[p].begin()); it;) {
				++sc;
				const Range<SeedLoc*> query_hits = *it.r;
				const Letter* seq = query_seqs.data(*query_hits.begin());
				if (!seed_is_complex(seq, shape, cut)) {
					++msc;
					qc += it.r->size();
					tc += it.s->size();
					for (SeedLoc* i = query_hits.begin(); i < query_hits.end(); ++i) {
						Letter* p = query_seqs.data(*i);
						*p |= SEED_MASK;
					}
					it.erase();
				}
				else
					++it;
			}
		}
		seed_count += sc;
		masked_seed_count += msc;
		query_count += qc;
		target_count += tc;
	};

	vector<std::thread> threads;
	for (int i = 0; i < config.threads_; ++i)
		threads.emplace_back(worker);
	for (auto& i : threads)
		i.join();
	timer.finish();
	verbose_stream << "Masked seeds: " << Util::String::ratio_percentage(masked_seed_count, seed_count) << endl;
	verbose_stream << "Masked positions (query): " << Util::String::ratio_percentage(query_count, query_seqs.letters()) << endl;
	verbose_stream << "Masked positions (target): " << Util::String::ratio_percentage(target_count, cfg.target->seqs().letters()) << endl;
}
