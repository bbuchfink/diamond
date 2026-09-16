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

#include "basic/config.h"
#include "search.h"
#include "util/math/integer.h"
#include "util/math/math.h"
#include "masking/def.h"
#include "basic/shape_config.h"
#include "align/def.h"
#include "align/extend.h"
#include "output/output_format.h"
#include "stats/cbs.h"

using std::vector;
using std::endl;
using std::map;
using std::string;
using std::prev;
using std::max;
using std::runtime_error;

namespace Search {

const double SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE = 0.15;

const map<Sensitivity, SensitivityTraits> sensitivity_traits = {
	{{ Sensitivity::FASTER, {
		true,      // support_query_indexed
		true,      // motif_masking
		50.0,      // freq_sd
		11,        // min_identities
		0,         // ungapped_evalue
		0,         // ungapped_evalue_short
		0,         // gapped_filter_evalue
		4,         // index_chunks
		16,        // query_bins
		nullptr,   // contiguous_seed
		0.9,       // seed_cut
		2.0,       // default_block_size
		murphy10,  // reduction
		0,         // minimizer_window
		21,        // sketch_size
		4,         // keyword_length
		12         // keyword_threshold
	}},
	{ Sensitivity::FAST, {
		true,      // support_query_indexed
		true,      // motif_masking
		50.0,      // freq_sd
		11,        // min_identities
		0,         // ungapped_evalue
		0,         // ungapped_evalue_short
		0,         // gapped_filter_evalue
		4,         // index_chunks
		16,        // query_bins
		nullptr,   // contiguous_seed
		0.9,       // seed_cut
		2.0,       // default_block_size
		murphy10,  // reduction
		0,         // minimizer_window
		0,         // sketch_size
		4,         // keyword_length
		12         // keyword_threshold
	}},
	{ Sensitivity::SHAPES6x10, {
		true,      // support_query_indexed
		true,      // motif_masking
		50.0,      // freq_sd
		11,        // min_identities
		0,         // ungapped_evalue
		0,         // ungapped_evalue_short
		0,         // gapped_filter_evalue
		4,         // index_chunks
		16,        // query_bins
		nullptr,   // contiguous_seed
		0.9,       // seed_cut
		2.0,       // default_block_size
		murphy10,  // reduction
		0,         // minimizer_window
		0,         // sketch_size
		4,         // keyword_length
		12         // keyword_threshold
	}},
	{ Sensitivity::SHAPES30x10, {
		true,      // support_query_indexed
		true,      // motif_masking
		50.0,      // freq_sd
		11,        // min_identities
		0,         // ungapped_evalue
		0,         // ungapped_evalue_short
		0,         // gapped_filter_evalue
		4,         // index_chunks
		16,        // query_bins
		nullptr,   // contiguous_seed
		0.9,       // seed_cut
		2.0,       // default_block_size
		murphy10,  // reduction
		0,         // minimizer_window
		0,         // sketch_size
		4,         // keyword_length
		12         // keyword_threshold
	}},
	{ Sensitivity::DEFAULT, {
		true,      // support_query_indexed
		true,      // motif_masking
		50.0,      // freq_sd
		11,        // min_identities
		10000,     // ungapped_evalue
		10000,     // ungapped_evalue_short
		0,         // gapped_filter_evalue
		4,         // index_chunks
		16,        // query_bins
		"111111",  // contiguous_seed
		0.8,       // seed_cut
		2.0,       // default_block_size
		murphy10,  // reduction
		0,         // minimizer_window
		0,         // sketch_size
		3,         // keyword_length
		8.4        // keyword_threshold
	}},
	{ Sensitivity::LINCLUST_40, {
		true,      // support_query_indexed
		true,      // motif_masking
		50.0,      // freq_sd
		11,        // min_identities
		0,         // ungapped_evalue
		0,         // ungapped_evalue_short
		0,         // gapped_filter_evalue
		4,         // index_chunks
		16,        // query_bins
		nullptr,   // contiguous_seed
		0.9,       // seed_cut
		2.0,       // default_block_size
		murphy10,  // reduction
		0,         // minimizer_window
		0,         // sketch_size
		4,         // keyword_length
		12         // keyword_threshold
	}},
	{ Sensitivity::LINCLUST_20, {
		true,      // support_query_indexed
		true,      // motif_masking
		50.0,      // freq_sd
		11,        // min_identities
		0,         // ungapped_evalue
		0,         // ungapped_evalue_short
		0,         // gapped_filter_evalue
		4,         // index_chunks
		16,        // query_bins
		nullptr,   // contiguous_seed
		0.9,       // seed_cut
		2.0,       // default_block_size
		murphy10,  // reduction
		0,         // minimizer_window
		0,         // sketch_size
		4,         // keyword_length
		12         // keyword_threshold
	}},
	{ Sensitivity::MID_SENSITIVE, {
		true,      // support_query_indexed
		true,      // motif_masking
		20.0,      // freq_sd
		11,        // min_identities
		10000,     // ungapped_evalue
		10000,     // ungapped_evalue_short
		0,         // gapped_filter_evalue
		4,         // index_chunks
		16,        // query_bins
		nullptr,   // contiguous_seed
		1.0,       // seed_cut
		2.0,       // default_block_size
		murphy10,  // reduction
		0,         // minimizer_window
		0,          // sketch_size
		3,         // keyword_length
		8.4        // keyword_threshold
	}},
	{ Sensitivity::SENSITIVE, {
		true,      // support_query_indexed
		true,      // motif_masking
		20.0,      // freq_sd
		11,        // min_identities
		10000,     // ungapped_evalue
		10000,     // ungapped_evalue_short
		1,         // gapped_filter_evalue
		4,         // index_chunks
		16,        // query_bins
		"11111",   // contiguous_seed
		1.0,       // seed_cut
		2.0,       // default_block_size
		murphy10,  // reduction
		0,         // minimizer_window
		0,          // sketch_size
		3,         // keyword_length
		8.4        // keyword_threshold
	}},
	{ Sensitivity::MORE_SENSITIVE, {
		true,      // support_query_indexed
		false,     // motif_masking
		200.0,     // freq_sd
		11,        // min_identities
		10000,     // ungapped_evalue
		10000,     // ungapped_evalue_short
		1,         // gapped_filter_evalue
		4,         // index_chunks
		16,        // query_bins
		"11111",   // contiguous_seed
		1.0,       // seed_cut
		2.0,       // default_block_size
		murphy10,  // reduction
		0,         // minimizer_window
		0,         // sketch_size
		3,         // keyword_length
		8.4        // keyword_threshold
	}},
	{ Sensitivity::VERY_SENSITIVE, {
		true,      // support_query_indexed
		false,     // motif_masking
		15.0,      // freq_sd
		9,         // min_identities
		100000,    // ungapped_evalue
		30000,     // ungapped_evalue_short
		1,         // gapped_filter_evalue
		1,         // index_chunks
		16,        // query_bins
		nullptr,   // contiguous_seed
		1.0,       // seed_cut
		0.4,       // default_block_size
		diamond9,  // reduction
		0,         // minimizer_window
		0,         // sketch_size
		3,         // keyword_length
		8.4        // keyword_threshold
	}},
	{ Sensitivity::ULTRA_SENSITIVE, {
		true,      // support_query_indexed
		false,     // motif_masking
		20.0,      // freq_sd
		9,         // min_identities
		300000,    // ungapped_evalue
		30000,     // ungapped_evalue_short
		1,         // gapped_filter_evalue
		1,         // index_chunks
		64,        // query_bins
		nullptr,   // contiguous_seed
		1.0,       // seed_cut
		0.4,       // default_block_size
		murphy10,  // reduction
		0,         // minimizer_window
		0,         // sketch_size
		3,         // keyword_length
		8.4        // keyword_threshold
	}},
} };

const map<Sensitivity, vector<Round>> iterated_sens{
	{ Sensitivity::FASTER,            { }},
	{ Sensitivity::FAST,            { {Sensitivity::FAST, true} }},
	{ Sensitivity::DEFAULT,         { {Sensitivity::FAST, true}, { Sensitivity::LINCLUST_40, true }}},
	{ Sensitivity::LINCLUST_40,     { {Sensitivity::FAST, true}, {Sensitivity::LINCLUST_40, true} }},
	{ Sensitivity::LINCLUST_20,     { {Sensitivity::FAST, true}, {Sensitivity::LINCLUST_20, true} }},
	{ Sensitivity::SHAPES30x10,     { {Sensitivity::FAST, true}, {Sensitivity::SHAPES30x10, true} }},
	{ Sensitivity::MID_SENSITIVE,   { {Sensitivity::FAST, true}, { Sensitivity::LINCLUST_40, true }, Sensitivity::DEFAULT}},
	{ Sensitivity::SENSITIVE,       { {Sensitivity::FAST, true}, { Sensitivity::LINCLUST_20, true }, Sensitivity::DEFAULT}},
	{ Sensitivity::MORE_SENSITIVE,  { {Sensitivity::FAST, true}, { Sensitivity::LINCLUST_20, true }, Sensitivity::DEFAULT}},
	{ Sensitivity::VERY_SENSITIVE,  { {Sensitivity::FAST, true}, { Sensitivity::LINCLUST_20, true }, Sensitivity::DEFAULT, Sensitivity::MORE_SENSITIVE}},
	{ Sensitivity::ULTRA_SENSITIVE, { {Sensitivity::FAST, true}, { Sensitivity::LINCLUST_20, true }, Sensitivity::DEFAULT, Sensitivity::MORE_SENSITIVE}}
};

const map<double, unsigned> approx_id_to_hamming_id {
	{ 50.0, 20 },
	{ 90.0, 30 }
};

static unsigned hamming_id_cutoff(double approx_id) {
	const auto it = approx_id_to_hamming_id.upper_bound(approx_id);
	return it == approx_id_to_hamming_id.begin() ? 0 : prev(it)->second;
}

const map<Sensitivity, vector<string>> shape_codes ={

	{{Sensitivity::DEFAULT, {
		"111101110111",
		"111011010010111" }},	// 2x10 iedera
	{Sensitivity::SENSITIVE, {
		"1011110111",
		"110100100010111",
		"11001011111",
		"101110001111",
		"11011101100001",
		"1111010010101",
		"111001001001011",
		"10101001101011",
		"111101010011",
		"1111000010000111",
		"1100011011011",
		"1101010000011011",
		"1110001010101001",
		"110011000110011",
		"11011010001101",
		"1101001100010011" }}, // 16x8 iedera
	{Sensitivity::MORE_SENSITIVE, {
		"1011110111",
		"110100100010111",
		"11001011111",
		"101110001111",
		"11011101100001",
		"1111010010101",
		"111001001001011",
		"10101001101011",
		"111101010011",
		"1111000010000111",
		"1100011011011",
		"1101010000011011",
		"1110001010101001",
		"110011000110011",
		"11011010001101",
		"1101001100010011" }}, // 16x8 iedera
	{ Sensitivity::VERY_SENSITIVE, {
		"11110111","11100100111","110010101011","11010001001011","10101100001101","110100100100011",
"1010010100010011","1100101000001011","11100000100010101","11000100010010011","11010000001000111",
"110001001000010011","1010001000100001011","1100010100000010011","1100100000101000011",
"1101000010000001011"
	}}, // 14x7
	{ Sensitivity::ULTRA_SENSITIVE, {
		"1111111",
		"11101111",
		"110011111",
		"110110111",
		"111111001",
		"1010111011",
		"1011110101",
		"1111000111",
		"10011110011",
		"10101101101",
		"10111010101",
		"11001010111",
		"11001100111",
		"11010101101",
		"11110001011",
		"100111010011",
		"101100110101",
		"101110000111",
		"110100101011",
		"110110001101",
		"111000110011",
		"1010001011011",
		"1010101000111",
		"1010110100011",
		"1100100110011",
		"1100101001011",
		"1101001100101",
		"1101010101001",
		"1110001010101",
		"1110010010011",
		"10100001101101",
		"11000100010111",
		"11010000100111",
		"11010100110001",
		"11101000011001",
		"11110000001101",
		"11110100000011",
		"101001000001111",
		"110000100101011",
		"110010010000111",
		"110101100001001",
		"110110000010011",
		"111001000100011",
		"111100000100101",
		"1000110010010101",
		"1001000100101101",
		"1001000110011001",
		"1010001001001011",
		"1010001010010011",
		"1010010001010101",
		"1010010100010011",
		"1010010101001001",
		"1010100000101011",
		"1010100011000101",
		"1011000010001011",
		"1100010000111001",
		"1100010010001011",
		"1100100001001011",
		"1100100100100011",
		"1100110000001101",
		"1101000100010011",
		"1101000110000101",
		"1110000001010011",
		"1110100000010101", // 64x7
}},
	{ Sensitivity::MID_SENSITIVE, {
		"11110110111",
		"1101100111101",
		"1110010101111",
		"11010101100111",
		"11101110001011",
		"1110100100010111",
		"1101000011010111",
		"1110011000011011"
}}, // 8x9
	{ Sensitivity::FAST, 
		{ "1101110101101111" } },
	{ Sensitivity::FASTER,
		{ "1101110101101111" } },
	{ Sensitivity::SHAPES6x10, {
"10111111111",
"111110110111",
"1101110111011",
"111111101011",
"1111011110011",
"111111100100011" } },
	{ Sensitivity::SHAPES30x10, {
		"10111111111",
		"111110110111",
		"1101110111011",
		"111111101011",
		"1111011110011",
		"111111100100011",
		"110111010011011",
"1111100110010011",
"11101100111101",
"111011011010101",
"11011010101111",
"11111110000010011",
"11011001100110011",
"101011100011111",
"111011111101",
"111110101100101",
"1111010101001011",
"11100111011001001",
"1110110001111001",
"110111011000010011",
"11001100101100111",
"11111000000111101",
"11011110011010001",
"110101101010011001",
"111010111000010101",
"1111101000100010011",
"11010100100111011",
"101001111100111",
"101110010001010111",
"11001101001011011"
	} },
	{ Sensitivity::LINCLUST_20, {
		"111111111111",
"1111111011111",
"1111110111111",
"11111111010111",
"11011101111111",
"11111011110111",
"11110011111111",
"11101111101111",
"11110111111011",
"110111110110111",
"111101111011011",
"111101100111111",
"111010111110111",
"111101011101111",
"111110110011111",
"111011101011111",
"111111010011111",
"111111001111011",
"111110101101111",
"111011110101111",
"1110101110011111",
"1111100110110111",
"1110111001101111",
"1111110010101111",
"1111001010111111",
"1110101101110111",
"1110110111001111",
"1110110101110111",
"1111010101101111",
"1111011011010111" }
},
{ Sensitivity::LINCLUST_40, {
		"111111111111",
"1111111011111",
"1111110111111",
"11111111010111",
"11011101111111",
"11111011110111",
"11110011111111",
"11101111101111",
"11110111111011",
"110111110110111",
"111101111011011",
"111101100111111",
"111010111110111",
"111101011101111",
"111110110011111" }
}
}
};

// Translated frames that repeat letter for letter are searched once. Frames are hashed, the first
// of each group is the representative, and the rest are skipped during seed enumeration; their hits
// are replayed from the representative in stage 2.
void build_frame_dups(Search::Config& cfg) {
	TaskTimer timer("Building frame duplicate map");
	const SequenceSet& seqs = cfg.query->seqs();
	const BlockId n = seqs.size();
	const size_t table_bits = std::max<size_t>(16, bit_length((int64_t)n) + 1), table_size = (size_t)1 << table_bits;
	std::vector<uint32_t> table(table_size, UINT32_MAX);
	std::vector<uint32_t> rep(n, UINT32_MAX);
	cfg.frame_skip.reset(new std::vector<bool>(n, false));
	std::vector<uint32_t> dup_count(n, 0);
	int64_t dups = 0, dup_letters = 0;

	for (BlockId i = 0; i < n; ++i) {
		const Sequence seq = seqs[i];
		uint64_t h = 1469598103934665603ull;
		for (Loc j = 0; j < seq.length(); ++j) {
			h ^= (unsigned char)letter_mask(seq[j]);
			h *= 1099511628211ull;
		}
		size_t slot = (size_t)(h >> (64 - table_bits));
		while (true) {
			const uint32_t e = table[slot];
			if (e == UINT32_MAX) {
				table[slot] = i;
				break;
			}
			const Sequence other = seqs[e];
			if (other.length() == seq.length() && std::equal(seq.data(), seq.end(), other.data())) {
				rep[i] = e;
				(*cfg.frame_skip)[i] = true;
				++dup_count[e];
				++dups;
				dup_letters += seq.length();
				break;
			}
			slot = (slot + 1) & (table_size - 1);
		}
	}

	cfg.frame_dup_begin.assign(n + 1, 0);
	for (BlockId i = 0; i < n; ++i)
		cfg.frame_dup_begin[i + 1] = cfg.frame_dup_begin[i] + dup_count[i];
	cfg.frame_dup_ids.resize(cfg.frame_dup_begin[n]);
	std::vector<uint32_t> pos(cfg.frame_dup_begin.begin(), cfg.frame_dup_begin.end() - 1);
	for (BlockId i = 0; i < n; ++i)
		if (rep[i] != UINT32_MAX)
			cfg.frame_dup_ids[pos[rep[i]]++] = i;

	timer.finish();
	*log_stream << "Duplicate frames = " << dups << '/' << n << " (" << (n ? 100.0 * dups / n : 0.0) << " %), letters = "
		<< dup_letters << '/' << seqs.letters() << std::endl;
}

int seedp_bits(int shape_weight, int threads, int index_chunks) {
	return max(max(bit_length(power((int64_t)Reduction::get_reduction().size(), (int64_t)shape_weight) - 1) - (int)sizeof(SeedOffset) * 8,
		bit_length((int64_t)threads * 4 * index_chunks - 1)), 8);
}

bool use_single_indexed(double coverage, size_t query_letters, size_t ref_letters)
{
	if (coverage >= SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE)
		return false;
	if (config.sensitivity >= Sensitivity::SENSITIVE) {
		return query_letters < 300000llu && query_letters * 20000llu < ref_letters;
	}
	else
		return query_letters < 3000000llu && query_letters * 2000llu < ref_letters;
}

MaskingAlgo soft_masking_algo(const SensitivityTraits& traits) {
	if (config.motif_masking.empty())
		return (!config.swipe_all && !config.freq_masking && traits.motif_masking) ? MaskingAlgo::MOTIF : MaskingAlgo::NONE;
	else {
		if (config.motif_masking == "0")
			return MaskingAlgo::NONE;
		else if (config.motif_masking == "1") {
			if (config.swipe_all)
				throw runtime_error("Soft masking is not supported for --swipe.");
			return MaskingAlgo::MOTIF;
		}
		else
			throw runtime_error("Permitted values for --motif-masking: 0, 1");
	}
}

void setup_search(Sensitivity sens, Search::Config& cfg)
{
	const SensitivityTraits& traits = sensitivity_traits.at(sens);
	config.sensitivity = sens;
	::Config::set_option(cfg.freq_sd, config.freq_sd_, 0.0, traits.freq_sd);
	::Config::set_option(cfg.hamming_filter_id, config.min_identities_, 0u, max(traits.min_identities, hamming_id_cutoff(config.approx_min_id.get(0.0))));
	::Config::set_option(cfg.ungapped_evalue, config.ungapped_evalue_, -1.0, traits.ungapped_evalue);
	::Config::set_option(cfg.ungapped_evalue_short, config.ungapped_evalue_short_, -1.0, traits.ungapped_evalue_short);
	::Config::set_option(cfg.gapped_filter_evalue, config.gapped_filter_evalue_, -1.0, traits.gapped_filter_evalue);
	if (config.query_bins.present())
		cfg.query_bins = (unsigned)config.query_bins.get_present();
	else
		cfg.query_bins = std::max((unsigned)std::round((double)config.threads_ / 8), traits.query_bins);
	::Config::set_option(cfg.minimizer_window, config.minimizer_window_, 0, traits.minimizer_window);
	::Config::set_option(cfg.sketch_size, config.sketch_size, 0, traits.sketch_size);

	cfg.keyword_length = traits.keyword_length;
	::Config::set_option(cfg.keyword_threshold, config.word_threshold, 0.0, traits.keyword_threshold);

	if (config.algo == ::Config::Algo::CTG_SEED) {
		if (!traits.contiguous_seed)
			throw runtime_error("Contiguous seed mode is not supported for this sensitivity setting.");
		if (sens == Sensitivity::DEFAULT)
			Reduction::set_reduction("KR EQ D N C G H F Y IV LM W P S T A");
		::shapes = ShapeConfig({ traits.contiguous_seed }, 0);
		Reduction::set_reduction(traits.reduction);
	}
	else {
		// The reduction has to be set before the shapes are built: Shape caches the letter
		// mask of the hashed seed encoding, whose field width is that of the reduction.
		Reduction::set_reduction(traits.reduction);
		::shapes = ShapeConfig(config.shape_mask.empty() ? shape_codes.at(sens) : config.shape_mask, config.shapes);
	}

	/* The on-disk seed index of --target-indexed is built over all seed positions
	   (HashedSeedSet), so subsampling the query would silently drop hits. */
	if (config.target_indexed) {
		cfg.minimizer_window = 0;
		cfg.sketch_size = 0;
	}
	if ((cfg.lin_stage1_target || config.lin_stage1_query) && shapes[0].weight_ < 10)
		throw runtime_error("Linearization is only supported for seed shapes of weight >= 10.");

	config.gapped_filter_diag_score = score_matrix.rawscore(config.gapped_filter_diag_bit_score);
	const double seed_cut = config.seed_cut_ == 0.0 ? traits.seed_cut : config.seed_cut_;
	cfg.seed_complexity_cut = seed_cut * Math::log(2.0) * ::shapes[0].weight_;
	cfg.soft_masking = soft_masking_algo(traits);
	if (!config.soft_masking.empty())
		cfg.soft_masking |= from_string<MaskingAlgo>(config.soft_masking);
	cfg.cutoff_table = Util::Scores::CutoffTable { cfg.ungapped_evalue };
	cfg.cutoff_table_short = Util::Scores::CutoffTable { cfg.ungapped_evalue_short };

	if (config.extension_mode.empty()) {
		if (config.global_ranking_targets || config.swipe_all)
			cfg.extension_mode = Extension::Mode::FULL;
		else
			cfg.extension_mode = Extension::default_ext_mode.at(sens);
	}
	else {
		cfg.extension_mode = from_string<Extension::Mode>(config.extension_mode);
		if (cfg.extension_mode != Extension::Mode::FULL) {
			if (config.global_ranking_targets)
				throw runtime_error("Global ranking only supports full matrix extension.");
			if (config.swipe_all)
				throw runtime_error("--swipe only supports full matrix extension.");
		}
	}

	if (cfg.extension_mode == Extension::Mode::FULL) {
		if (config.frame_shift > 0)
			throw runtime_error("Frameshift alignment does not support full matrix extension.");
	}

	// The new extension pipeline computes the extensions using the anchored swipe, which
	// computes no traceback and uses the unadjusted scoring matrix.
	if (config.new_extension_pipeline) {
		if (config.comp_based_stats_.get(Stats::DEFAULT_CBS) != Stats::CBS::DISABLED)
			throw runtime_error("The new extension pipeline requires --comp-based-stats 0.");
		if (!flag_only(cfg.output_format->hsp_values, HspValues::COORDS))
			throw runtime_error("The new extension pipeline only supports output fields that do not require a traceback.");
		if (cfg.extension_mode == Extension::Mode::FULL)
			throw runtime_error("The new extension pipeline does not support full matrix extension.");
	}

	if(!config.aln_out.empty() && config.parallel_tmpdir.empty())
		throw runtime_error("Alignment output to file is not supported without --parallel-tmpdir.");
}

}