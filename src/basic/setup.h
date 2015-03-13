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

#include <omp.h>
#include <boost/iostreams/tee.hpp>
#include "options.h"
#include "../util/system.h"

using std::cout;
using std::endl;

void setup(const string &command, int ac, const char **av)
{
	namespace io = boost::iostreams;
	namespace po = program_options;
	if(po::debug_log) {
		io::tee_filter<io::file_sink> t (io::file_sink("diamond.log", std::ios_base::out | std::ios_base::app));
		verbose_stream.push(t);
		log_stream.push(t);
		log_stream.push(cout);
	} else
		log_stream.push(io::null_sink ());
	if(po::verbose || po::debug_log)
		verbose_stream.push(cout);
	else
		verbose_stream.push(io::null_sink ());

	log_stream << "Command line: ";
	for(int i=0;i<ac;++i)
		log_stream << av[i] << ' ';
	log_stream << endl;

	verbose_stream << Const::program_name << " v" << Const::version_string << "." << Const::build_version << endl;
#ifndef NDEBUG
	verbose_stream << "Assertions enabled." << endl;
#endif
	po::set_option(po::threads_, boost::thread::hardware_concurrency());
	omp_set_num_threads(po::threads_);
	omp_set_dynamic(0);
	verbose_stream << "#Threads = " << omp_get_max_threads() << endl;

	if(command == "makedb")
		po::command = po::makedb;
	else if(command == "blastx")
		po::command = po::blastx;
	else if(command == "blastp")
		po::command = po::blastp;
	else if(command == "blastn")
		po::command = po::blastn;
	else if(command == "view")
		po::command = po::view;
	else
		po::command = po::invalid;

	if(sequence_type() == amino_acid) {
		if(po::gap_open == -1)
			po::gap_open = 11;
		if(po::gap_extend == -1)
			po::gap_extend = 1;
		score_matrix::instance = auto_ptr<score_matrix> (new score_matrix(po::matrix,
				po::gap_open,
				po::gap_extend,
				po::reward,
				po::penalty,
				Amino_acid ()));
		score_matrix::get().print<Amino_acid>();
	} else {
#ifdef EXTRA
		if(po::gap_open == -1)
			po::gap_open = 5;
		if(po::gap_extend == -1)
			po::gap_extend = 2;
		score_matrix::instance = auto_ptr<score_matrix> (new score_matrix(po::matrix,
				po::gap_open,
				po::gap_extend,
				po::reward,
				po::penalty,
				Nucleotide ()));
		score_matrix::get().print<Nucleotide>();
#endif
	}
	verbose_stream << "Gap open penalty = " << po::gap_open << endl;
	verbose_stream << "Gap extension penalty = " << po::gap_extend << endl;

	if(po::seg == "" && po::command == po::blastx)
		po::seg = "yes";
	verbose_stream << "Seg masking = " << (po::seg == "yes") << endl;

	po::have_ssse3 = check_SSSE3();
	if(po::have_ssse3)
		verbose_stream << "SSSE3 enabled." << endl;
	if(po::debug_log) {
		copy_file(log_stream, "/etc/issue");
		copy_file(log_stream, "/proc/cpuinfo");
		copy_file(log_stream, "/proc/meminfo");
	}
}

template<typename _val>
void setup_search_params(pair<size_t,size_t> query_len_bounds, size_t chunk_db_letters)
{
	namespace po = program_options;
	if(po::aligner_mode == po::sensitive) {
		po::set_option(po::hit_cap, 256u);
	} else if (po::aligner_mode == po::fast) {
		po::set_option(po::hit_cap, 32u);
	}

	const double b = po::min_bit_score == 0 ? score_matrix::get().bitscore(po::max_evalue, ref_header.letters, query_len_bounds.first) : po::min_bit_score;

	if(query_len_bounds.second <= 40) {
		po::set_option(po::min_identities, 10u);
		po::set_option(po::min_ungapped_raw_score, score_matrix::get().rawscore(std::min(27.0, b)));
	} else {
		po::set_option(po::min_identities, 9u);
		po::set_option(po::min_ungapped_raw_score, score_matrix::get().rawscore(std::min(23.0, b)));
	}

	if(query_len_bounds.second <= 80) {
		const int band = po::read_padding<_val>(query_len_bounds.second);
		po::set_option(po::window, (unsigned)(query_len_bounds.second + band));
		po::set_option(po::hit_band, band);
		po::set_option(po::min_hit_score, score_matrix::get().rawscore(b));
	} else {
		po::set_option(po::window, 40u);
		po::set_option(po::hit_band, 5);
		po::set_option(po::min_hit_score, score_matrix::get().rawscore(std::min(29.0, b)));
	}
	log_stream << "Query len bounds " << query_len_bounds.first << ' ' << query_len_bounds.second << endl;
	log_stream << "Search parameters " << po::min_ungapped_raw_score << ' ' << po::min_hit_score << ' ' << po::hit_cap << endl;
}

template<>
void setup_search_params<Amino_acid>(pair<size_t,size_t> query_len_bounds, size_t chunk_db_letters)
{
	namespace po = program_options;
	if(po::aligner_mode == po::sensitive) {
		po::set_option(po::hit_cap, std::max(256u, (unsigned)(chunk_db_letters/8735437)));
	} else if (po::aligner_mode == po::fast) {
		po::set_option(po::hit_cap, std::max(128u, (unsigned)(chunk_db_letters/17470874)));
	}

	const double b = po::min_bit_score == 0 ? score_matrix::get().bitscore(po::max_evalue, ref_header.letters, query_len_bounds.first) : po::min_bit_score;

	if(query_len_bounds.second <= 40) {
		po::set_option(po::min_identities, 10u);
		po::set_option(po::min_ungapped_raw_score, score_matrix::get().rawscore(std::min(27.0, b)));
	} else {
		po::set_option(po::min_identities, 9u);
		po::set_option(po::min_ungapped_raw_score, score_matrix::get().rawscore(std::min(23.0, b)));
	}

	if(query_len_bounds.second <= 80) {
		const int band = po::read_padding<Amino_acid>(query_len_bounds.second);
		po::set_option(po::window, (unsigned)(query_len_bounds.second + band));
		po::set_option(po::hit_band, band);
		po::set_option(po::min_hit_score, score_matrix::get().rawscore(b));
	} else {
		po::set_option(po::window, 40u);
		po::set_option(po::hit_band, 5);
		po::set_option(po::min_hit_score, score_matrix::get().rawscore(std::min(29.0, b)));
	}
	log_stream << "Query len bounds " << query_len_bounds.first << ' ' << query_len_bounds.second << endl;
	log_stream << "Search parameters " << po::min_ungapped_raw_score << ' ' << po::min_hit_score << ' ' << po::hit_cap << endl;
}

#endif /* SETUP_H_ */
