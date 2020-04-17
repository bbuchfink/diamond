#include <string>
#include <vector>
#include <stdexcept>
#include <queue>
#include <utility>
#include <iostream>
#include "../basic/config.h"
#include "../util/io/text_input_file.h"
#include "../util/seq_file_format.h"
#include "../util/util.h"
#include "../util/log_stream.h"
#include "../util/string/tokenizer.h"

using std::string;
using std::vector;
using std::endl;

struct Record {

	enum { BLANK = -1 };

	Record():
		query_id(BLANK)
	{
	}

	Record(TextInputFile &file, size_t file_id) :
		file(file_id)
	{
		file.getline();
		if (file.line.empty()) {
			query_id = BLANK;
			return;
		}
		Util::String::Tokenizer(file.line, "\t") >> query_id >> query_acc >> subject_acc >> evalue;
	}

	bool operator<(const Record &r) const {
		return query_id > r.query_id || (query_id == r.query_id && (evalue > r.evalue || (evalue == r.evalue && file > r.file)));
	}
	
	bool blank() const {
		return query_id == BLANK;
	}

	friend std::ostream& operator<<(std::ostream& os, const Record& r) {
		os << r.query_acc << '\t' << r.subject_acc << '\t' << r.evalue << endl;
		return os;
	}

	int query_id;
	string query_acc, subject_acc;
	double evalue;
	size_t file;

};

void merge_tsv() {
	if (config.input_ref_file.empty())
		throw std::runtime_error("Missing parameter --in");
	task_timer timer("Processing input");
	
	const size_t n = config.input_ref_file.size();
	message_stream << "#Input files: " << n << endl;
	vector<TextInputFile> files;
	std::priority_queue<Record> queue;
	size_t records = 0;
	for (size_t i = 0; i < n; ++i) {
		files.emplace_back(config.input_ref_file[i]);
		queue.emplace(files.back(), i);
		++records;
	}

	while (!queue.empty()) {
		if (queue.top().blank()) {
			queue.pop();
			continue;
		}
		std::cout << queue.top() << endl;
		const size_t file = queue.top().file;

		Record r;
		do
			r = Record(files[file], file);
		while (!r.blank() && queue.top().query_acc == r.query_acc && queue.top().subject_acc == r.subject_acc);
		queue.pop();
		if (!r.blank()) {
			queue.push(std::move(r));
			++records;
		}
	}

	for (auto &f : files)
		f.close();
	message_stream << "#Records: " << records << endl;
}