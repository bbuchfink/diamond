#include <unordered_map>
#include <string>
#include <vector>
#include <stdexcept>
#include <queue>
#include <utility>
#include <iostream>
#include "tsv_record.h"
#include "../basic/config.h"
#include "../util/io/text_input_file.h"
#include "../util/seq_file_format.h"
#include "../util/util.h"
#include "../util/log_stream.h"

using std::string;
using std::vector;
using std::endl;

static std::unordered_map<string, int> acc2idx;

struct Record : public TSVRecord {

	Record(TSVRecord &&r, size_t file) :
		TSVRecord(r),
		file(file)
	{
		auto it = acc2idx.find(qseqid);
		if (it == acc2idx.end())
			throw std::runtime_error("Query accession not mapped.");
		query_id = it->second;
	}

	bool operator<(const Record &r) const {
		return query_id > r.query_id || (query_id == r.query_id && (evalue > r.evalue || (evalue == r.evalue && bitscore < r.bitscore)));
	}

	int query_id;
	size_t file;

};

void merge_tsv() {
	task_timer timer("Loading query file");
	TextInputFile query_file(config.query_file);
	
	FASTA_format format;
	vector<Letter> seq;
	string id;
	while (format.get_seq(id, seq, query_file))
		acc2idx[blast_id(id)] = acc2idx.size();
	message_stream << "#Queries: " << acc2idx.size() << endl;

	timer.go("Processing input");
	if (config.input_ref_file.empty())
		throw std::runtime_error("Missing parameter --in");
	const size_t n = config.input_ref_file.size();
	message_stream << "#Input files: " << n << endl;
	vector<TextInputFile> files;
	TSVRecord r;
	std::priority_queue<Record> queue;
	size_t records = 0;
	for (size_t i = 0; i < n; ++i) {
		files.emplace_back(config.input_ref_file[i]);
		files.back().getline();
		files.back() >> r;
		if (!r.blank()) {
			queue.emplace(std::move(r), i);
			++records;
		}
	}

	while (!queue.empty()) {
		std::cout << queue.top() << endl;
		const size_t file = queue.top().file;
		do
			files[file] >> r;
		while (!r.blank() && queue.top().qseqid == r.qseqid && queue.top().sseqid == r.sseqid);
		queue.pop();
		if (!r.blank()) {
			queue.emplace(std::move(r), file);
			++records;
		}
	}

	for (auto &f : files)
		f.close();
	message_stream << "#Records: " << records << endl;
}