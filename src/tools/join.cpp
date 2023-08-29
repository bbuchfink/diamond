#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include "../util/io/text_input_file.h"
#include "../basic/config.h"
#include "../util/tsv/tsv.h"
#include "../util/string/tokenizer.h"
#include "../util/log_stream.h"

using std::cout;
using std::endl;
using std::string;
using namespace Util::Tsv;

void join() {
	TextInputFile file1(config.file1);
	TextInputFile file2(config.file2);
	std::unordered_map<string, string> values;
	std::unordered_set<string> keys;
	file1.getline();
	const size_t key2 = config.key2 - 1;

	if (column_count(file1.line) == 1) {

		do {
			keys.insert(file1.line);
		} while (file1.getline(), !file1.line.empty() || !file1.eof());

		message_stream << "#Keys: " << keys.size() << endl;

		while (file2.getline(), !file2.line.empty() || !file2.eof()) {
			if (keys.find(column(file2.line, key2)) != keys.end())
				cout << file2.line << endl;
		}

	}
	else {

		do {
			values[Util::Tsv::column(file1.line, 0)] = Util::Tsv::column(file1.line, 1);
		} while (file1.getline(), !file1.line.empty() || !file1.eof());

		string key;
		while (file2.getline(), !file2.line.empty() || !file2.eof()) {
			Util::String::Tokenizer tok(file2.line, "\t");
			tok >> key;
			cout << key << '\t' << values.at(key) << '\t' << tok.ptr() << std::endl;
		}

	}

	file1.close();
	file2.close();
}