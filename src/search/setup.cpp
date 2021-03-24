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
		"10111111",
"111101011",
"1101100111",
"10110011011",
"11000111011",
"110100101011",
"110110000111",
"10101001101001",
"11001010010011",
"11010000110011",
"11010011000101",
"11100100010101",
"101010100010011",
"1110000001001011",

		/*"11111011",
		"101101111",
		"1101010111",
		"1110100111",
		"11001101011",
		"11100011011",
		"1110010001011",
		"11010001000111",
		"110100001010011",
		"1010100100001011",
		"1101000001000111",
		"1110001000010101",
		"1100010010010011",
		"1100101000001011"*/
	}}, // 16x7 SpEED 0.35 35
	{ Sensitivity::ULTRA_SENSITIVE, {
		"1111111",
"10111111",
"11011111",
"110111101",
"111101011",
"111101101",
"1010011111",
"1011110011",
"1101100111",
"1110111001",
"10110011011",
"10111010101",
"11000111011",
"11011101001",
"11110100011",
"101001001111",
"110000111101",
"110011001011",
"110011100011",
"110100101011",
"110110000111",
"1010001010111",
"1011010000111",
"1101001001011",
"10000011100111",
"10011000100111",
"10101001101001",
"10101010011001",
"11001010010011",
"11010000110011",
"11010011000101",
"11011100000101",
"11100100010101",
"100010010110101",
"100011010011001",
"100110001001101",
"101000100110101",
"101010001100101",
"101010100010011",
"101100100001011",
"110001100100101",
"110010000010111",
"110100010001101",
"110100100100011",
"111001000001101",
"1000011100010101",
"1001000110100011",
"1001001000101011",
"1001010100110001",
"1010000101010101",
"1010010101000011",
"1010100010110001",
"1010101000010011",
"1010101110000001",
"1100000110100101",
"1100010110000011",
"1100100100011001",
"1100110000000111",
"1101101000001001",
"1110000000110101",
"1110000001001011",
"1110000010010101",
"1110001010001001",
"1111000000010011"

		/*"1111111",
		"11100100111",
		"101110001101",
		"110101001011",
		"1110000101011",
		"1101100010101",
		"1110000110101",
		"1110000100111",
		"11010001100101",
		"10100100110011",
		"11010010001011",
		"10101000100111",
		"11000101010011",
		"11000110001011",
		"11001010100011",
		"11010000011101",
		"110010010010101",
		"110010001001011",
		"110001001100011",
		"110100010001011",
		"110100100100101",
		"101100010000111",
		"101001010001011",
		"110010010010011",
		"101011000010011",
		"110010100001011",
		"101100000110011",
		"110001001001101",
		"101010100100011",
		"110100010100011",
		"111000100001101",
		"110110000001101",
		"1100010000101011",
		"1100001010000111",
		"1001010010000111",
		"1100101000001011",
		"1010011000001011",
		"1110000001010011",
		"1010100010100011",
		"1010010011000011",
		"1101001000100011",
		"1010001100000111",
		"1010010000110011",
		"1101001000010101",
		"1100100001010011",
		"1110001000100101",
		"1100010100010011",
		"1101000100100011",
		"1101000010010011",
		"1011000001100011",
		"1010100001001101",
		"1100100101000011",
		"1101000011001001",
		"1100010100100101",
		"1100100100010101",
		"1100001100001011",
		"1011001000010011",
		"1001000101001011",
		"1010100010000111",
		"1100010001010101",
		"1100100001100101",
		"1010100100001011",
		"1100101000101001",
		"1100010010001101" // 64x7*/
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