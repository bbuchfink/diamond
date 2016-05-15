/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef SETUP_H_
#define SETUP_H_

#include "../basic/config.h"

/*void setup_search_params(pair<size_t,size_t> query_len_bounds, size_t chunk_db_letters)
{
	namespace po = program_options;
	if(po::aligner_mode == po::sensitive) {
		po::set_option(po::hit_cap, 256u);
	} else if (po::aligner_mode == po::fast) {
		po::set_option(po::hit_cap, 32u);
	}

	const double b = po::min_bit_score == 0 ? score_matrix::get().bitscore(po::max_evalue, ref_header.letters, query_len_bounds.first) : po::min_bit_score;

	po::set_option(po::min_identities, 18u);
	po::min_ungapped_raw_score = score_matrix::get().rawscore(std::min(po::min_ungapped_raw_score == 0 ? 19.0 : po::min_ungapped_raw_score, b));

	po::set_option(po::window, 40u);
	po::set_option(po::hit_band, 5);
	po::min_hit_score = score_matrix::get().rawscore(std::min(po::min_hit_score == 0 ? 19.0 : po::min_hit_score, b));

	log_stream << "Query len bounds " << query_len_bounds.first << ' ' << query_len_bounds.second << endl;
	log_stream << "Minimum bit score = " << b << endl;
	log_stream << "Search parameters " << po::min_ungapped_raw_score << ' ' << po::min_hit_score << ' ' << po::hit_cap << endl;
}*/

void setup_search_params(pair<size_t,size_t> query_len_bounds, size_t chunk_db_letters)
{
	if(config.mode_sensitive) {
		Config::set_option(config.hit_cap, std::max(256u, (unsigned)(chunk_db_letters/8735437)));
	} else {
		Config::set_option(config.hit_cap, std::max(128u, (unsigned)(chunk_db_letters/17470874)));
	}

	const double b = config.min_bit_score == 0 ? score_matrix.bitscore(config.max_evalue, ref_header.letters, (unsigned)query_len_bounds.first) : config.min_bit_score;

	if(query_len_bounds.second <= 40) {
		Config::set_option(config.min_identities, 10u);
		Config::set_option(config.min_ungapped_raw_score, score_matrix.rawscore(std::min(27.0, b)));
	} else {
		Config::set_option(config.min_identities, 9u);
		Config::set_option(config.min_ungapped_raw_score, score_matrix.rawscore(std::min(23.0, b)));
	}

	if(query_len_bounds.second <= 80) {
		const int band = config.read_padding(query_len_bounds.second);
		Config::set_option(config.window, (unsigned)(query_len_bounds.second + band));
		Config::set_option(config.hit_band, band);
		Config::set_option(config.min_hit_score, score_matrix.rawscore(b));
	} else {
		Config::set_option(config.window, 40u);
		Config::set_option(config.hit_band, 5);
		Config::set_option(config.min_hit_score, score_matrix.rawscore(std::min(29.0, b)));
	}
	log_stream << "Query len bounds " << query_len_bounds.first << ' ' << query_len_bounds.second << endl;
	log_stream << "Search parameters " << config.min_ungapped_raw_score << ' ' << config.min_hit_score << ' ' << config.hit_cap << endl;
}

#endif /* SETUP_H_ */
