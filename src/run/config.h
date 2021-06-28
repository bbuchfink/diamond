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
#include "../util/data_structures/bit_vector.h"
#include "../util/scores/cutoff_table.h"

struct SequenceFile;
struct Consumer;
struct TextInputFile;
struct Block;
struct TaxonomyNodes;
enum class Sensitivity;
enum class SeedEncoding;
template<typename T> struct AsyncBuffer;

struct Async;
template<typename T, size_t E, typename Sync> struct Deque;

namespace Extension { namespace GlobalRanking {
	struct Hit;
}}

namespace Search {

struct Hit;

struct Config {

	using RankingTable = std::vector<Extension::GlobalRanking::Hit>;
	using RankingBuffer = Deque<Search::Hit, 28, Async>;

	Config();
	void free();
	~Config();

	bool                                       self;
	std::vector<Sensitivity>                   sensitivity;
	SeedEncoding                               seed_encoding;
	bool                                       track_aligned_queries;
	double                                     freq_sd;
	unsigned                                   hamming_filter_id;
	double                                     ungapped_evalue;
	double                                     ungapped_evalue_short;
	double                                     gapped_filter_evalue;
	unsigned                                   index_chunks;
	unsigned                                   query_bins;

	std::shared_ptr<SequenceFile>              db;
	std::shared_ptr<std::list<TextInputFile>>  query_file;
	std::shared_ptr<Consumer>                  out;
	std::shared_ptr<BitVector>                 db_filter;
	TaxonomyNodes*                             taxon_nodes;
	std::vector<std::string>*                  taxonomy_scientific_names;

	std::unique_ptr<Block>                     query, target;
	std::unique_ptr<std::vector<bool>>         query_skip;
	std::unique_ptr<AsyncBuffer<Hit>>          seed_hit_buf;
	std::unique_ptr<RankingBuffer>             global_ranking_buffer;
	std::unique_ptr<RankingTable>              ranking_table;

	uint64_t db_seqs, db_letters, ref_blocks;
	Util::Scores::CutoffTable cutoff_gapped1, cutoff_gapped2;
	Util::Scores::CutoffTable2D cutoff_gapped1_new, cutoff_gapped2_new;

	size_t                                     iteration_query_aligned;

	bool iterated() const {
		return sensitivity.size() > 1;
	}

};

}
