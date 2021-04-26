/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
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

#include "../data/reference.h"
#include "../basic/config.h"
#include "seed_complexity.h"
#include "search.h"

using std::endl;
using std::map;

double SeedComplexity::prob_[AMINO_ACID_COUNT];
const double SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE = 0.15;

const map<Sensitivity, SensitivityTraits> sensitivity_traits {
	//                               qidx   freqsd minid ug_ev   ug_ev_s gf_ev  idx_chunk qbins
	{ Sensitivity::FAST,            {true,  50.0,  11,   10000,  10000,  0,     4,        16 }},
	{ Sensitivity::DEFAULT,         {true,  50.0,  11,   10000,  10000,  0,     4,        16 }},
	{ Sensitivity::MID_SENSITIVE,   {true,  20.0,  11,   10000,  10000,  0,     4,        16 }},
	{ Sensitivity::SENSITIVE,       {true,  20.0,  11,   10000,  10000,  1,     4,        16 }},
	{ Sensitivity::MORE_SENSITIVE,  {true,  200.0, 11,   10000,  10000,  1,     4,        16 }},
	{ Sensitivity::VERY_SENSITIVE,  {true,  15.0,  9,    100000, 30000,  1,     1,        16 }},
	{ Sensitivity::ULTRA_SENSITIVE, {true,  20.0,  9,    300000, 30000,  1,     1,        64 }}
};

static const map<Sensitivity, vector<string>> shape_codes = {

	{Sensitivity::DEFAULT, {
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
		{ "1101110101101111" } }

};

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

void setup_search()
{
	const Sensitivity sens = config.sensitivity;
	const SensitivityTraits& traits = sensitivity_traits.at(sens);
	Config::set_option(config.freq_sd, traits.freq_sd);
	Config::set_option(config.min_identities, traits.min_identities);
	Config::set_option(config.ungapped_evalue, traits.ungapped_evalue, -1.0);
	Config::set_option(config.ungapped_evalue_short, traits.ungapped_evalue_short, -1.0);
	Config::set_option(config.gapped_filter_evalue, traits.gapped_filter_evalue, -1.0);
	Config::set_option(config.lowmem, traits.index_chunks);
	Config::set_option(config.query_bins, traits.query_bins);
	
	if (config.algo == Config::Algo::QUERY_INDEXED) {
		config.lowmem = 1;
		if(!traits.support_query_indexed)
			throw std::runtime_error("Query-indexed mode is not supported for this sensitivity setting.");
	}
	
	::shapes = ShapeConfig(config.shape_mask.empty() ? shape_codes.at(sens) : config.shape_mask, config.shapes);

	print_warnings();

	if (config.command != Config::blastp && config.command != Config::blastx)
		return;

	SeedComplexity::init(Reduction::reduction);
	config.gapped_filter_diag_score = score_matrix.rawscore(config.gapped_filter_diag_bit_score);

	verbose_stream << "Seed frequency SD: " << config.freq_sd << endl;
	verbose_stream << "Shape configuration: " << ::shapes << endl;
}