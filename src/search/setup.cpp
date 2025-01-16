/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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

#include "basic/config.h"
#include "search.h"
#include "util/math/integer.h"
#include "masking/def.h"
#include "basic/shape_config.h"

using std::vector;
using std::endl;
using std::map;
using std::string;
using std::prev;
using std::max;

namespace Search {

const double SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE = 0.15;

Reduction murphy10("A KR EDNQ C G H ILVM FYW P ST");
Reduction steinegger12("AST C DN EQ FY G H IV KR LM P W");
Reduction dna("A C G T");

const map<Sensitivity, SensitivityTraits> sensitivity_traits = {
	//                               qidx   motifm freqsd minid ug_ev   ug_ev_s gf_ev  idx_chunk qbins ctg_seed  seed_cut block_size reduction min_window sketch
	{{ Sensitivity::FASTER,         {true,  true,  50.0,  11,   0,      0,      0,     4,        16,   nullptr,  0.9,     2.0,       murphy10, 0,         21 }},
	{ Sensitivity::FAST,            {true,  true,  50.0,  11,   0,      0,      0,     4,        16,   nullptr,  0.9,     2.0,       murphy10, 0,         0 }},
	{ Sensitivity::SHAPES30x10,     {true,  true,  50.0,  11,   0,      0,      0,     4,        16,   nullptr,  0.9,     2.0,       murphy10, 0,         0 }},
	{ Sensitivity::DEFAULT,         {true,  true,  50.0,  11,   10000,  10000,  0,     4,        16,   "111111", 0.8,     2.0,       murphy10, 0,         0 }},
	{ Sensitivity::LINCLUST_20,     {true,  true,  50.0,  11,   0,      0,      0,     4,        16,   nullptr,  0.9,     2.0,       murphy10, 0,         0 }},
	{ Sensitivity::MID_SENSITIVE,   {true,  true,  20.0,  11,   10000,  10000,  0,     4,        16,   nullptr,  1.0,     2.0,       murphy10, 0,         0 }},
	{ Sensitivity::SENSITIVE,       {true,  true,  20.0,  11,   10000,  10000,  1,     4,        16,   "11111",  1.0,     2.0,       murphy10, 0,         0 }},
	{ Sensitivity::MORE_SENSITIVE,  {true,  false, 200.0, 11,   10000,  10000,  1,     4,        16,   "11111",  1.0,     2.0,       murphy10, 0,         0 }},
	{ Sensitivity::VERY_SENSITIVE,  {true,  false, 15.0,  9,    100000, 30000,  1,     1,        16,   nullptr,  1.0,     0.4,       murphy10, 0,         0 }},
	{ Sensitivity::ULTRA_SENSITIVE, {true,  false, 20.0,  9,    300000, 30000,  1,     1,        64,   nullptr,  1.0,     0.4,       murphy10, 0,         0 }},
} };

const map<Sensitivity, vector<Round>> iterated_sens{
	{ Sensitivity::FASTER,            { }},
	{ Sensitivity::FAST,            { {Sensitivity::FAST, true} }},
	{ Sensitivity::DEFAULT,         { {Sensitivity::FAST, true}, Sensitivity::FAST }},
	{ Sensitivity::LINCLUST_20,     { {Sensitivity::FAST, true}, {Sensitivity::LINCLUST_20, true} }},
	{ Sensitivity::SHAPES30x10,     { {Sensitivity::FAST, true}, {Sensitivity::SHAPES30x10, true} }},
	{ Sensitivity::MID_SENSITIVE,   { {Sensitivity::FAST, true}, Sensitivity::FAST, Sensitivity::DEFAULT}},
	{ Sensitivity::SENSITIVE,       { {Sensitivity::FAST, true}, Sensitivity::FAST, Sensitivity::DEFAULT}},
	{ Sensitivity::MORE_SENSITIVE,  { {Sensitivity::FAST, true}, Sensitivity::FAST, Sensitivity::DEFAULT}},
	{ Sensitivity::VERY_SENSITIVE,  { {Sensitivity::FAST, true}, Sensitivity::FAST, Sensitivity::DEFAULT, Sensitivity::MORE_SENSITIVE}},
	{ Sensitivity::ULTRA_SENSITIVE, { {Sensitivity::FAST, true}, Sensitivity::FAST, Sensitivity::DEFAULT, Sensitivity::MORE_SENSITIVE}}
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
	{Sensitivity::VERY_SENSITIVE, {
		"11101111",
		"110110111",
		"111111001",
		"1010111011",
		"11110001011",
		"110100101011",
		"110110001101",
		"1010101000111",
		"1100101001011",
		"1101010101001",
		"1110010010011",
		"110110000010011",
		"111001000100011",
		"1101000100010011",
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
}}
};

int seedp_bits(int shape_weight, int threads, int index_chunks) {
	return max(max(bit_length(power((int64_t)Reduction::reduction.size(), (int64_t)shape_weight) - 1) - (int)sizeof(SeedOffset) * 8,
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

bool keep_target_id(const Search::Config& cfg) {
#ifdef HIT_KEEP_TARGET_ID
	return true;
#else
	return cfg.min_length_ratio != 0.0 || config.global_ranking_targets;
#endif
}

MaskingAlgo soft_masking_algo(const SensitivityTraits& traits) {
	if (config.motif_masking.empty())
		return (!config.swipe_all && !config.freq_masking && traits.motif_masking) ? MaskingAlgo::MOTIF : MaskingAlgo::NONE;
	else {
		if (config.motif_masking == "0")
			return MaskingAlgo::NONE;
		else if (config.motif_masking == "1") {
			if (config.swipe_all)
				throw std::runtime_error("Soft masking is not supported for --swipe.");
			return MaskingAlgo::MOTIF;
		}
		else
			throw std::runtime_error("Permitted values for --motif-masking: 0, 1");
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
	::Config::set_option(cfg.query_bins, config.query_bins_, 0u, traits.query_bins);
	::Config::set_option(cfg.minimizer_window, config.minimizer_window_, 0, traits.minimizer_window);
	::Config::set_option(cfg.sketch_size, config.sketch_size, 0, traits.sketch_size);
	
	if (config.algo == ::Config::Algo::CTG_SEED) {
		if (!traits.contiguous_seed)
			throw std::runtime_error("Contiguous seed mode is not supported for this sensitivity setting.");
		if (sens == Sensitivity::DEFAULT)
			Reduction::reduction = Reduction("KR EQ D N C G H F Y IV LM W P S T A");
		::shapes = ShapeConfig({ traits.contiguous_seed }, 0);
	}
	else
		::shapes = ShapeConfig(config.shape_mask.empty() ? shape_codes.at(sens) : config.shape_mask, config.shapes);

	config.gapped_filter_diag_score = score_matrix.rawscore(config.gapped_filter_diag_bit_score);
	const double seed_cut = config.seed_cut_ == 0.0 ? traits.seed_cut : config.seed_cut_;
	cfg.seed_complexity_cut = seed_cut * std::log(2.0) * ::shapes[0].weight_;
	cfg.soft_masking = soft_masking_algo(traits);
	if (!config.soft_masking.empty())
		cfg.soft_masking |= from_string<MaskingAlgo>(config.soft_masking);
	cfg.cutoff_table = { cfg.ungapped_evalue };
	cfg.cutoff_table_short = { cfg.ungapped_evalue_short };
}

}
