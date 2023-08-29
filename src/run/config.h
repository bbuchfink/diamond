/****
DIAMOND protein aligner
Copyright (C) 2016-2021 Max Planck Society for the Advancement of Science e.V.
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

#pragma once
#include <list>
#include <memory>
#include <mutex>
#include "../util/data_structures/bit_vector.h"
#include "../util/scores/cutoff_table.h"
#include "../stats/dna_scoring/build_score.h"

struct SequenceFile;
struct Consumer;
struct TextInputFile;
struct Block;
struct TaxonomyNodes;
struct ThreadPool;
struct OutputFormat;
enum class Sensitivity;
enum class SeedEncoding;
enum class MaskingAlgo;
template<typename T> struct AsyncBuffer;

struct Async;
template<typename T, size_t E, typename Sync> struct Deque;

namespace Extension {
	enum class Mode;
	namespace GlobalRanking {
	struct Hit;
}}
namespace Dna{
    class Index;
}
namespace Search {

struct Round {
	Round(Sensitivity sens, bool lin = false):
		sensitivity(sens),
		linearize(lin)
	{}
	bool operator<(Round r) const {
		return sensitivity < r.sensitivity || (sensitivity == r.sensitivity && linearize && !r.linearize);
	}
	bool operator==(Round r) const {
		return sensitivity == r.sensitivity && linearize == r.linearize;
	}
	Sensitivity sensitivity;
	bool linearize;
};

struct Hit;

struct Config {

	using RankingTable = std::vector<Extension::GlobalRanking::Hit>;
	using RankingBuffer = Deque<Search::Hit, 28, Async>;

	Config();
	void free();
	~Config();

	bool                                       self;
	std::vector<Round>                         sensitivity;
	SeedEncoding                               seed_encoding;
	MaskingAlgo                                query_masking;
	MaskingAlgo                                target_masking;
	MaskingAlgo                                soft_masking;
	Extension::Mode                            extension_mode;
	double                                     seed_complexity_cut;
	bool                                       lazy_masking;
	bool                                       track_aligned_queries;
	double                                     freq_sd;
	Loc                                        minimizer_window;
	bool                                       lin_stage1_target;
	unsigned                                   hamming_filter_id;
	double                                     ungapped_evalue;
	double                                     ungapped_evalue_short;
	double                                     gapped_filter_evalue;
	unsigned                                   index_chunks;
	unsigned                                   query_bins;
	int64_t                                    max_target_seqs;
	std::unique_ptr<OutputFormat>              output_format;

	std::shared_ptr<SequenceFile>              db;
	std::shared_ptr<SequenceFile>              query_file;
	std::shared_ptr<Consumer>                  out;
	std::shared_ptr<BitVector>                 db_filter;

	std::shared_ptr<Block>                     query, target;
	std::unique_ptr<std::vector<bool>>         query_skip;
	std::unique_ptr<AsyncBuffer<Hit>>          seed_hit_buf;
	std::unique_ptr<RankingBuffer>             global_ranking_buffer;
	std::unique_ptr<RankingTable>              ranking_table;
    std::unique_ptr<Stats::Blastn_Score>       score_builder;
#ifdef WITH_DNA
    std::unique_ptr<Dna::Index>                dna_ref_index;
#endif

    int                                        current_query_block;
	int                                        current_ref_block;
	bool                                       blocked_processing;
	std::vector<bool>                          aligned_targets;
	std::mutex                                 aligned_targets_mtx;

    uint64_t db_seqs, db_letters, ref_blocks;
#ifdef UNGAPPED_SPOUGE
	const Util::Scores::CutoffTable2D cutoff_table;
#else
	Util::Scores::CutoffTable cutoff_table, cutoff_table_short;
#endif
	Util::Scores::CutoffTable cutoff_gapped1, cutoff_gapped2;
	Util::Scores::CutoffTable2D cutoff_gapped1_new, cutoff_gapped2_new;

	BlockId                                    iteration_query_aligned;

	std::unique_ptr<ThreadPool>                thread_pool;

	bool iterated() const {
		return sensitivity.size() > 1;
	}

};

}
