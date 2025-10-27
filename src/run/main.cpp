/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#include <iostream>
#include "basic/config.h"
#include "tools.h"
#include "workflow.h"
#include "cluster/cluster_registry.h"
#ifdef WITH_MCL
#include "contrib/mcl/recursive_parser.h"
#endif
#include "data/dmnd/dmnd.h"
#include "util/command_line_parser.h"
#include "util/log_stream.h"
#include "util/system/system.h"

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
//void roc();
//void roc_id();
void makeindex();
void find_shapes();
void hash_seqs();
void list_seeds();
void greedy_vertex_cover();
void merge_daa();
#ifdef EXTRA
void word_count();
void cut();
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
			Search::run();
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
			Workflow::Cluster::ClusterRegistry::get(config.cluster_algo.get("cascaded"))->run();
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
		case Config::find_shapes:
			find_shapes();
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
			greedy_vertex_cover();
			break;
		case Config::CLUSTER_REASSIGN:
			Cluster::reassign();
			break;
		case Config::RECLUSTER:
			Cluster::recluster();
			break;
		case Config::MERGE_DAA:
			merge_daa();
			break;
#ifdef EXTRA
        case Config::blastn:
            Search::run();
            break;
		case Config::WORD_COUNT:
			word_count();
			break;
		case Config::CUT:
			cut();
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
		cerr << "Failed to allocate sufficient memory. Please refer to the online wiki for instructions on memory usage." << endl;
		log_stream << "Error: " << e.what() << endl;
		return 1;
	}
	catch (const FileOpenException&) {
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