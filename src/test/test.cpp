/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

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
#include <string>
#include <limits.h>
#include <string.h>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <list>
#include "../util/io/temp_file.h"
#include "../util/io/text_input_file.h"
#include "test.h"
#include "../util/sequence/sequence.h"
#include "../util/log_stream.h"
#include "../data/reference.h"
#include "../basic/statistics.h"
#include "../run/workflow.h"
#include "../util/util.h"
#include "../util/string/string.h"
#include "../util/system/system.h"
#include "../data/dmnd/dmnd.h"

using std::endl;
using std::string;
using std::vector;
using std::cout;
using std::list;
using std::shared_ptr;

namespace Test {

size_t run_testcase(size_t i, shared_ptr<DatabaseFile> &db, shared_ptr<list<TextInputFile>>& query_file, size_t max_width, bool bootstrap, bool log, bool to_cout) {
	vector<string> args = tokenize(test_cases[i].command_line, " ");
	args.emplace(args.begin(), "diamond");
	if (log)
		args.push_back("--log");
	config = Config((int)args.size(), charp_array(args.begin(), args.end()).data(), false);
	statistics.reset();
	query_file->front().rewind();

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

int run() {
	const bool bootstrap = config.bootstrap, log = config.debug_log, to_cout = config.output_file == "stdout";
	task_timer timer("Generating test dataset");
	TempFile proteins;
	for (size_t i = 0; i < seqs.size(); ++i)
		Util::Seq::format(Sequence::from_string(seqs[i].second.c_str()), seqs[i].first.c_str(), nullptr, proteins, "fasta", amino_acid_traits);
	shared_ptr<list<TextInputFile>> query_file(new list<TextInputFile>);
	query_file->emplace_back(proteins);
	timer.finish();

	config.command = Config::makedb;
	TempFile *db_file;
	DatabaseFile::make_db(&db_file, query_file.get());
	shared_ptr<DatabaseFile> db(new DatabaseFile(*db_file));

	const size_t n = test_cases.size(),
		max_width = std::accumulate(test_cases.begin(), test_cases.end(), (size_t)0, [](size_t l, const TestCase& t) { return std::max(l, strlen(t.desc)); });
	size_t passed = 0;
	for (size_t i = 0; i < n; ++i)
		passed += run_testcase(i, db, query_file, max_width, bootstrap, log, to_cout);

	cout << endl << "#Test cases passed: " << passed << '/' << n << endl; // << endl;
	
	query_file->front().close_and_delete();
	db->close();
	delete db_file;
	return passed == n ? 0 : 1;
}

}