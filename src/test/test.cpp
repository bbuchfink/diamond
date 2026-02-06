/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#include <iostream>
#include <string>
#include <string.h>
#include <numeric>
#include <iomanip>
#include <list>
#include "util/io/temp_file.h"
#include "test.h"
#include "util/log_stream.h"
#include "basic/statistics.h"
#include "run/workflow.h"
#include "util/string/string.h"
#include "util/system/system.h"
#include "basic/config.h"
#include "data/fasta/fasta_file.h"
#include "util/command_line_parser.h"

using std::endl;
using std::string;
using std::vector;
using std::cout;
using std::list;
using std::shared_ptr;

int run_queue_stress_test();
void filestack();

namespace Test {

static size_t run_testcase(size_t i, shared_ptr<SequenceFile> &db, shared_ptr<SequenceFile>& query_file, size_t max_width, bool bootstrap, bool log, bool to_cout) {
	vector<string> args = Util::String::tokenize(test_cases[i].command_line, " ");
	args.emplace(args.begin(), "diamond");
	if (log)
		args.push_back("--log");
	CommandLineParser parser;
	config = Config((int)args.size(), charp_array(args.begin(), args.end()).data(), false, parser);
	statistics.reset();
	query_file->set_seqinfo_ptr(0);
	db->set_seqinfo_ptr(0);

	if (to_cout) {
		Search::run(db, query_file);
		return 0;
	}
	
	shared_ptr<TempFile> output_file(new TempFile(!bootstrap));

	Search::run(db, query_file, output_file);

	InputFile out_in(*output_file);
	uint64_t hash = out_in.hash();

	if (bootstrap)
		out_in.close();
	else
		out_in.close_and_delete();

	if (bootstrap)
		cout << "0x" << std::hex << hash << ',' << endl;
	else {
		const bool passed = hash == ref_hashes[i];
		cout << std::setw(max_width) << std::left << test_cases[i].desc << " [ ";
		set_color(passed ? Color::GREEN : Color::RED);
		cout << (passed ? "Passed" : "Failed");
		reset_color();
		cout << " ]" << endl;
		return passed ? 1 : 0;
	}
	return 0;
}

static void load_seqs(SequenceFile& file) {
	file.init_write();
	for (size_t i = 0; i < seqs.size(); ++i)
		file.write_seq(Sequence::from_string(seqs[i].second.c_str()), seqs[i].first.c_str());
}

int run() {
	const bool bootstrap = config.bootstrap, log = config.debug_log, to_cout = config.output_file == "stdout";
	//filestack();
	run_queue_stress_test();
	TaskTimer timer("Generating test dataset");
	FastaFile proteins("test1", true, FastaFile::WriteAccess());
	
	shared_ptr<SequenceFile> query_file(new FastaFile("test2", true, FastaFile::WriteAccess(), SequenceFile::Flags::ALL)), db(new FastaFile("test3", true, FastaFile::WriteAccess(), SequenceFile::Flags::ALL));
	load_seqs(*query_file);
	load_seqs(*db);
	timer.finish();

	const size_t n = test_cases.size(),
		max_width = std::accumulate(test_cases.begin(), test_cases.end(), (size_t)0, [](size_t l, const TestCase& t) { return std::max(l, strlen(t.desc)); });
	size_t passed = 0;
	for (size_t i = 0; i < n; ++i)
		passed += run_testcase(i, db, query_file, max_width, bootstrap, log, to_cout);

	cout << endl << "#Test cases passed: " << passed << '/' << n << endl;
	
	query_file->close();
	db->close();
	return passed == n ? 0 : 1;
}

}