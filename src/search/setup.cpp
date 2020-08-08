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

double SeedComplexity::prob_[AMINO_ACID_COUNT];
const double SINGLE_INDEXED_SEED_SPACE_MAX_COVERAGE = 0.15;

void setup_search_cont()
{
	if (config.sensitivity >= Sensitivity::VERY_SENSITIVE || config.sensitivity == Sensitivity::MID_SENSITIVE)
		return;
	unsigned index_mode;
	Reduction::reduction = Reduction("A KR EDNQ C G H ILVM FYW P ST");
	if (config.sensitivity == Sensitivity::SENSITIVE || config.sensitivity == Sensitivity::MORE_SENSITIVE)
		index_mode = 11;
	else {
		index_mode = 10;
		Reduction::reduction = Reduction("KR EQ D N C G H F Y IV LM W P S T A");
	}
	::shapes = shape_config(index_mode, 1, vector<string>());
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

void setup_search()
{
	if (config.sensitivity == Sensitivity::ULTRA_SENSITIVE) {
		Config::set_option(config.freq_sd, 20.0);
		Config::set_option(config.min_identities, 9u);
		Config::set_option(config.ungapped_evalue, 300000.0, -1.0);
		Config::set_option(config.gapped_filter_evalue, 10.0, -1.0);
		Config::set_option(config.lowmem, 1u);
		Config::set_option(config.query_bins, 64u);
	} else if (config.sensitivity == Sensitivity::VERY_SENSITIVE) {
		Config::set_option(config.freq_sd, 15.0);
		Config::set_option(config.min_identities, 9u);
		Config::set_option(config.ungapped_evalue, 100000.0, -1.0);
		Config::set_option(config.gapped_filter_evalue, 10.0, -1.0);
		Config::set_option(config.lowmem, 1u);
		Config::set_option(config.query_bins, 16u);
	} else if (config.sensitivity == Sensitivity::MORE_SENSITIVE) {
		Config::set_option(config.freq_sd, 200.0);
		Config::set_option(config.min_identities, 11u);
		Config::set_option(config.ungapped_evalue, 10000.0, -1.0);
		Config::set_option(config.gapped_filter_evalue, 10.0, -1.0);
		Config::set_option(config.lowmem, 4u);
		Config::set_option(config.query_bins, 16u);
	}
	else if (config.sensitivity == Sensitivity::SENSITIVE) {
		Config::set_option(config.freq_sd, 20.0);
		Config::set_option(config.min_identities, 11u);
		Config::set_option(config.ungapped_evalue, 10000.0, -1.0);
		Config::set_option(config.gapped_filter_evalue, 10.0, -1.0);
		Config::set_option(config.lowmem, 4u);
		Config::set_option(config.query_bins, 16u);
	}
	else if (config.sensitivity == Sensitivity::MID_SENSITIVE) {
		Config::set_option(config.freq_sd, 20.0);
		Config::set_option(config.min_identities, 11u);
		Config::set_option(config.ungapped_evalue, 10000.0, -1.0);
		Config::set_option(config.gapped_filter_evalue, 0.0, -1.0);
		Config::set_option(config.lowmem, 4u);
		Config::set_option(config.query_bins, 16u);
	}
	else {
		Config::set_option(config.freq_sd, 50.0);
		Config::set_option(config.min_identities, 11u);
		Config::set_option(config.ungapped_evalue, 10000.0, -1.0);
		Config::set_option(config.lowmem, 4u);
		Config::set_option(config.query_bins, 16u);
	}
	
	if(config.algo==Config::query_indexed)
		config.lowmem = 1;
	else {
		switch (config.sensitivity) {
		case Sensitivity::ULTRA_SENSITIVE:
			Config::set_option(config.index_mode, 13u);
			break;
		case Sensitivity::VERY_SENSITIVE:
			Config::set_option(config.index_mode, 12u);
			break;
		case Sensitivity::MORE_SENSITIVE:
		case Sensitivity::SENSITIVE:
			Config::set_option(config.index_mode, 9u);
			break;
		case Sensitivity::MID_SENSITIVE:
			Config::set_option(config.index_mode, 15u);
			break;
		case Sensitivity::FAST:
			Config::set_option(config.index_mode, 8u);
			break;
		}
		Reduction::reduction = Reduction("A KR EDNQ C G H ILVM FYW P ST");
		::shapes = shape_config(config.index_mode, config.shapes, config.shape_mask);
	}

	print_warnings();

	SeedComplexity::init(Reduction::reduction);

	message_stream << "Algorithm: " << (config.algo == Config::double_indexed ? "Double-indexed" : "Query-indexed") << endl;
	verbose_stream << "Reduction: " << Reduction::reduction << endl;

	verbose_stream << "Seed frequency SD: " << config.freq_sd << endl;
	verbose_stream << "Shape configuration: " << ::shapes << endl;
}