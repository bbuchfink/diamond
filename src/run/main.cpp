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

#include <iostream>
#include "basic/config.h"
#include "workflow.h"
#ifdef WITH_MCL
#include "contrib/mcl/recursive_parser.h"
#endif
#include "legacy/dmnd/dmnd.h"
#include "util/command_line_parser.h"
#include "util/log_stream.h"
#include "util/system/system.h"
#include "tools/tools.h"

using std::cout;
using std::cerr;
using std::endl;
using std::runtime_error;
using std::string;

void opt();
void run_masker();
void fastq2fasta();
void view_daa();
void db_info();
void benchmark_sw();
void db_annot_stats();
void info();
void seed_stat();
void pairwise();
void reverse();
void makeindex();
void hash_seqs();
void list_seeds();
void merge_daa();
void multinode();
#ifdef EXTRA
namespace Cluster {
#ifdef WITH_FAMSA
	void profile_recluster();
#endif
}
#endif

void split();
namespace Benchmark { void benchmark(); }
namespace Test { int run();
}
namespace Cluster {
void reassign();
void realign();
void recluster();
namespace Incremental {
}}

int main(int ac, const char* av[])
{
	std::unique_ptr<std::vector<BitVector>> target_seed_hits;
	std::set_terminate([]() noexcept {
		std::abort();
		// std::_Exit(EXIT_FAILURE);
		// std::quick_exit(EXIT_FAILURE);
		});
	try {
		init_motif_table();
		CommandLineParser parser;
		config = Config(ac, av, true, parser);
		switch (config.command) {
		case Config::help:
			break;
		case Config::version:
			cout << Const::program_name << " version " << Const::version_string << endl;
			break;
		case Config::makedb:
			DatabaseFile::make_db();
			break;
		case Config::blastp:
		case Config::blastx:
			Search::run(target_seed_hits);
			break;
		case Config::view:
			if (!config.daa_file.empty())
				view_daa();
			else
				throw std::runtime_error("The view command requires a DAA (option -a) input file.");
			break;
		case Config::getseq:
			get_seq();
			break;
		case Config::random_seqs:
			random_seqs();
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
		case Config::info:
			info();
			break;
		case Config::smith_waterman:
			pairwise();
			break;
		case Config::cluster:
		case Config::DEEPCLUST:
		case Config::LINCLUST:
#ifdef WITH_MCL
			// Why is cluster_similarity not set at the end of the Config constructor?
			if(!config.cluster_similarity.empty()){
				string expression = RecursiveParser::clean_expression(&config.cluster_similarity);
				RecursiveParser rp(nullptr, expression.c_str());
				try{
					rp.evaluate();
				}
				catch (const runtime_error& e){
					message_stream << "Could not evaluate the expression: " << config.cluster_similarity << endl;
					throw e;
				}
			}
#endif
			multinode();
			break;
		case Config::benchmark:
			Benchmark::benchmark();
			break;
		case Config::split:
			split();
			break;
		case Config::regression_test:
			return Test::run();
			break;
		case Config::reverse_seqs:
			reverse();
			break;
		case Config::roc:
			throw std::runtime_error("Deprecated command: roc");
			break;
		case Config::rocid:
			throw std::runtime_error("Deprecated command: rocid");
			break;
		case Config::makeidx:
			makeindex();
			break;
		case Config::HASH_SEQS:
			hash_seqs();
			break;
		case Config::prep_db:
			set_color(Color::YELLOW, true);
			cerr << "Warning: prepdb is deprecated since v2.1.14 and no longer needed to use BLAST databases. No action was taken." << endl;
			reset_color(true);
			break;
		case Config::LIST_SEEDS:
			list_seeds();
			break;
		case Config::CLUSTER_REALIGN:
			Cluster::realign();
			break;
		case Config::GREEDY_VERTEX_COVER:
			GVC::greedy_vertex_cover();
			break;
		case Config::CLUSTER_REASSIGN:
			set_color(Color::YELLOW, true);
			cerr << "Reassign has been temporarily removed for v2.2.1. No action was taken." << endl;
			reset_color(true);
			//Cluster::reassign();
			break;
		case Config::RECLUSTER:
			set_color(Color::YELLOW, true);
			cerr << "Recluster has been temporarily removed for v2.1.25. No action was taken." << endl;
			reset_color(true);
			//Cluster::recluster();
			break;
		case Config::MERGE_DAA:
			merge_daa();
			break;
#ifdef EXTRA
        case Config::blastn:
            Search::run(target_seed_hits);
            break;
#ifdef WITH_FAMSA
		case Config::PROFILE_RECLUSTER:
			Cluster::profile_recluster();
			break;
#endif
#endif
		default:
			return 1;
		}
	}
	catch (const std::bad_alloc &e) {
		cleanup();
		cerr << "Failed to allocate sufficient memory. Please refer to the online wiki for instructions on memory usage." << endl;
		*log_stream << "Error: " << e.what() << endl;
		return 1;
	}
	catch (const FileOpenException&) {
		cleanup();
		return 1;
	} catch(const std::exception& e) {
		cleanup();
        cerr << "Error: " << e.what() << endl;
		*log_stream << "Error: " << e.what() << endl;
        return 1;
    }
    catch(...) {
		cleanup();
        cerr << "Exception of unknown type!" << endl;
        return 1;
    }

	cleanup();
    return 0;
}