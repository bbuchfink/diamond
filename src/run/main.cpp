/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#include <iostream>
#include "../basic/config.h"
#include "tools.h"
#include "../extra/compare.h"

using std::cout;
using std::cerr;
using std::endl;

void master_thread_di();
void model_seqs();
void opt();
void run_masker();
void fastq2fasta();
void view();
void make_db();
void db_info();
void test_main();
void benchmark_sw();
void test_io();
void db_annot_stats();
void read_sim();
void info();
void seed_stat();
void pairwise();
void protein_snps();

int main(int ac, const char* av[])
{
	try {

		check_simd();
		config = Config(ac, av);

		switch (config.command) {
		case Config::help:
			break;
		case Config::version:
			cout << Const::program_name << " version " << Const::version_string << endl;
			break;
		case Config::makedb:
			make_db();
			break;
		case Config::blastp:
		case Config::blastx:
			master_thread_di();
			break;
		case Config::view:
			view();
			break;
		case Config::getseq:
			get_seq();
			break;
		case Config::random_seqs:
			random_seqs();
			break;
		case Config::compare:
			compare();
			break;
		case Config::sort:
			sort_file();
			break;
		case Config::roc:
			roc();
			break;
		case Config::db_stat:
			db_stat();
			break;
		case Config::model_sim:
			model_sim();
			break;
		case Config::match_file_stat:
			match_file_stat();
			break;
		case Config::model_seqs:
			model_seqs();
			break;
		case Config::opt:
			opt();
			break;
		case Config::mask:
			run_masker();
			break;
		case Config::fastq2fasta:
			fastq2fasta();
			break;
		case Config::dbinfo:
			db_info();
			break;
		case Config::test_io:
			test_io();
			break;
		case Config::db_annot_stats:
			db_annot_stats();
			break;
		case Config::read_sim:
			read_sim();
			break;
		case Config::info:
			info();
			break;
		case Config::smith_waterman:
			pairwise();
			break;
		case Config::protein_snps:
			protein_snps();
			break;
#ifdef EXTRA
		case Config::test_extra:
			test_main();
			break;
		case Config::benchmark:
			benchmark_sw();
			break;
		case Config::seed_stat:
			seed_stat();
			break;
#endif
		default:
			return 1;
		}
	}
	catch(const std::bad_alloc &e) {
		cerr << "Failed to allocate sufficient memory. Please refer to the manual for instructions on memory usage." << endl;
		log_stream << "Error: " << e.what() << endl;
		return 1;
	} catch(const std::exception& e) {
        cerr << "Error: " << e.what() << endl;
        log_stream << "Error: " << e.what() << endl;
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!" << endl;
        return 1;
    }

    return 0;
}
