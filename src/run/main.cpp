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

#include <iostream>
#include "../basic/config.h"
#include "../output/view.h"
#include "tools.h"
#include "../extra/compare.h"

using std::cout;
using std::cerr;
using std::endl;

void run_mapper();
void master_thread_di();
void model_seqs();
void opt();

int main(int ac, const char* av[])
{

	try {

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
			if (config.algo == Config::subject_indexed)
				run_mapper();
			else
				master_thread_di();
			break;
		case Config::view:
			view();
			break;
		case Config::getseq:
			get_seq();
			break;
		case Config::benchmark:
			benchmark_sw();
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
		default:
			return 1;
		}
	}
	catch(std::bad_alloc &e) {
		cerr << "Failed to allocate sufficient memory. Please refer to the manual for instructions on memory usage." << endl;
		log_stream << "Error: " << e.what() << endl;
		return 1;
	} catch(std::exception& e) {
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
