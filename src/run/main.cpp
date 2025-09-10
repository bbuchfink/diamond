/****
DIAMOND protein aligner
Copyright (C) 2013-2022 Max Planck Society for the Advancement of Science e.V.
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

#include <iostream>
#ifdef WITH_MIMALLOC
#include <mimalloc-2.0/mimalloc-new-delete.h>
#endif
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
void translate();
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
void length_sort();
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
		case Config::translate:
			translate();
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
		case Config::LENGTH_SORT:
			length_sort();
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
		cerr << "Failed to allocate sufficient memory. Please refer to the manual for instructions on memory usage." << endl;
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
