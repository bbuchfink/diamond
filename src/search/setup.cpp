/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#include "align_range.h"
#include "../data/reference.h"
#include "../basic/config.h"

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
