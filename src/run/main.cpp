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
#include "../basic/const.h"
#include "../util/log_stream.h"
#include "../data/reference.h"
#include "../run/make_db.h"
#include "../run/master_thread.h"
#include "../output/view.h"
#include "tools.h"

#ifdef EXTRA
#include "../../extra/test_sw.h"
#include "../../extra/test_io.h"
#include "../../extra/stat.h"
#endif

#ifdef ENABLE_STAT
#include "../stat/compare.h"
//#include "../stat/random_seq.h"
#endif

using std::cout;
using std::cerr;
using std::endl;

int main(int ac, const char* av[])
{

	try {

		config = Config(ac, av);

		switch(config.command) {
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
			master_thread();
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
		#ifdef EXTRA
        else if (command == "stat" && vm.count("match1"))
        	if(program_options::db_type == "nucl")
        		blast_stat<Nucleotide>(vm.count("tab") > 0);
        	else
        		blast_stat<Amino_acid>(vm.count("tab") > 0);
        //else if (command == "comp" && vm.count("query") && vm.count("match1") && vm.count("match2"))
        	//compare();
        //else if (command == "random" && vm.count("query") && vm.count("out"))
        	//random_seq();
		#endif
		#ifdef EXTRA
        else if (command == "test")
        	test_io();
		#endif
		default:
        	return 1;
	}
	}
	catch(std::bad_alloc &e) {
		cerr << "Failed to allocate sufficient memory. Please refer to the readme for instructions on memory usage." << endl;
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
