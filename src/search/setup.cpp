/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include "align_range.h"
#include "../data/reference.h"
#include "../basic/config.h"

void setup_search_cont()
{
	unsigned index_mode;
	if (config.mode_more_sensitive) {
		index_mode = 11;
	}
	else if (config.mode_sensitive) {
		index_mode = 11;
	}
	else {
		index_mode = 10;
		Reduction::reduction = Reduction("KR EQ D N C G H F Y IV LM W P S T A");
	}
	::shapes = shape_config(index_mode, 1, vector<string>());
}

void setup_search()
{
	if (config.algo == Config::double_indexed) {
		if (config.mode_more_sensitive) {
			Config::set_option(config.index_mode, 9u);
			Config::set_option(config.freq_sd, 200.0);
		}
		else if (config.mode_sensitive) {
			Config::set_option(config.index_mode, 9u);
			Config::set_option(config.freq_sd, 10.0);
		}
		else {
			Config::set_option(config.index_mode, 8u);
			Config::set_option(config.freq_sd, 50.0);
		}
		Reduction::reduction = Reduction("A KR EDNQ C G H ILVM FYW P ST");
		::shapes = shape_config(config.index_mode, config.shapes, config.shape_mask);
	}
	else {
		if (config.mode_more_sensitive) {
			Config::set_option(config.freq_sd, 200.0);
		}
		else if (config.mode_sensitive) {
			Config::set_option(config.freq_sd, 20.0);
		}
		else {
			Config::set_option(config.freq_sd, 50.0);
		}
		config.lowmem = 1;
	}

	message_stream << "Algorithm: " << (config.algo == Config::double_indexed ? "Double-indexed" : "Query-indexed") << endl;
	verbose_stream << "Reduction: " << Reduction::reduction << endl;

	verbose_stream << "Seed frequency SD: " << config.freq_sd << endl;
	verbose_stream << "Shape configuration: " << ::shapes << endl;
	config.seed_anchor = std::min(::shapes[0].length_ - 1, 8u);
}

void setup_search_params(pair<size_t, size_t> query_len_bounds, size_t chunk_db_letters)
{
	const double b = config.min_bit_score == 0 ? score_matrix.bitscore(config.max_evalue, ref_header.letters, (unsigned)query_len_bounds.first) : config.min_bit_score;

	if (config.mode_very_sensitive) {
		Reduction::reduction = Reduction("A KR EDNQ C G H ILVM FYW P ST"); // murphy.10
		Config::set_option(config.index_mode, 4u);
		::shapes = shape_config(config.index_mode, config.shapes, config.shape_mask);
		config.seed_anchor = std::min(::shapes[0].length_ - 1, 8u);
		Config::set_option(config.min_identities, 9u);
		Config::set_option(config.min_ungapped_score, 19.0);
		Config::set_option(config.window, 60u);
		Config::set_option(config.hit_band, 8);
		Config::set_option(config.min_hit_score, 23.0);
	}
	else {

		Config::set_option(config.min_identities, 11u);
		if (query_len_bounds.second <= 40) {
			//Config::set_option(config.min_identities, 10u);
			Config::set_option(config.min_ungapped_score, std::min(27.0, b));
		}
		else {
			//Config::set_option(config.min_identities, 9u);
			Config::set_option(config.min_ungapped_score, std::min(23.0, b));
		}

		if (query_len_bounds.second <= 80) {
			const int band = config.read_padding(query_len_bounds.second);
			Config::set_option(config.window, (unsigned)(query_len_bounds.second + band));
			Config::set_option(config.hit_band, band);
			Config::set_option(config.min_hit_score, b);
		}
		else {
			Config::set_option(config.window, 40u);
			Config::set_option(config.hit_band, 5);
			Config::set_option(config.min_hit_score, std::min(29.0, b));
		}
	}

	config.min_ungapped_raw_score = score_matrix.rawscore(config.min_ungapped_score);
	config.min_hit_raw_score = score_matrix.rawscore(config.min_hit_score);
	log_stream << "Query len bounds " << query_len_bounds.first << ' ' << query_len_bounds.second << endl;
	log_stream << "Search parameters " << config.min_ungapped_raw_score << ' ' << config.min_hit_score << ' ' << config.hit_cap << endl;
}
