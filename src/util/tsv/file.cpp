#include "file.h"
#include "../enum.h"
#include "construct.h"
#include "../data_structures/reorder_queue.h"

using std::vector;
using std::function;
using std::string;
using std::runtime_error;

namespace Util { namespace Tsv {

static Flags get_flags(Flags flags) {
	if (flag_any(flags, Flags::TEMP)) {
		if (flag_any(flags, Flags::WRITE))
			throw runtime_error("Write-only temp file.");
		flags |= Flags::READ_WRITE | Flags::OVERWRITE;
	}
	return flags;
}

File::File(const Schema& schema, const char* file_name, Flags flags, const Config& config) :
	flags_(get_flags(flags)),
	schema_(schema),
	config_(config),
	out_file_(flag_any(flags_, Flags::WRITE | Flags::READ_WRITE) ? (flag_any(flags_, Flags::TEMP) ? new TempFile : new OutputFile(file_name, Compressor::NONE, flag_any(flags_, Flags::OVERWRITE) ? "w+b" : "r+b")) : nullptr),
	file_(flag_any(flags_, Flags::READ_WRITE) ? (flag_any(flags_, Flags::TEMP) ? new TextInputFile(*(TempFile*)out_file_.get()) : new TextInputFile(*out_file_)) : (flag_any(flags_, Flags::WRITE) ? nullptr : new TextInputFile(file_name))),
	record_id_(0)
{
	if (flag_any(flags_, Flags::RECORD_ID_COLUMN) && schema.front() != Type::INT64)
		throw runtime_error("Schema does not contain record_id column.");
}

File::File(const Schema& schema, const std::string& file_name, Flags flags, const Config& config):
	File(schema, file_name.c_str(), flags, config)
{}

void File::rewind() {
	if (out_file_)
		out_file_->rewind();
	if (file_)
		file_->rewind();
	record_id_ = 0;
}

File::~File() {
	if (file_ && file_->temp_file)
		file_->close_and_delete();
	else {
		if (out_file_)
			out_file_->close();
		else
			file_->close();
	}
}

void File::read(int threads, std::function<void(int64_t chunk, const Table&)>& callback) {
	function<void(int64_t, const char*, const char*)> f([this, &callback](int64_t chunk, const char* begin, const char* end) {
		Table table(schema_);
		table.append(begin, end);
		callback(chunk, table);
		});
	read(threads, f);
}

Table File::read(int threads) {
	Table table(schema_);
	function<void(const char*, const char*)> f([this](const char* begin, const char* end) {
		Table table(this->schema());
		LineIterator it(begin, end);
		while (it.good()) {
			const string line = *it;
			table.push_back(line.begin(), line.end());
			++it;
		}
	});
	rewind();
	while(file_->getline(), (!file_->eof() || !file_->line.empty())) {
		table.push_back(file_->line.cbegin(), file_->line.cend());
	}
	return table;
}

Table File::read(int64_t max_records, int threads) {
	Table table(schema_);
	while (table.size() < max_records && (file_->getline(), (!file_->eof() || !file_->line.empty()))) {
		table.push_back(file_->line.cbegin(), file_->line.cend(), flag_any(flags_, Flags::RECORD_ID_COLUMN) ? record_id_ : -1);
		++record_id_;
	}
	return table;
}

void File::write_record(int i) {
	if (i != (flag_any(flags_, Flags::RECORD_ID_COLUMN) ? 1 : 0))
		throw runtime_error("write_record with insufficient field count.");
	write_buf_ << '\n';
	out_file_->write(write_buf_.data(), write_buf_.size());
	write_buf_.clear();
}

void File::write(const Record& record) {
	record.write(write_buf_);
	out_file_->write(write_buf_.data(), write_buf_.size());
	write_buf_.clear();
}

void File::write(const Table& table) {
	for (int64_t i = 0; i < table.size(); ++i)
		write(table[i]);
}

File* File::map(int threads, std::function<Table(const Record&)>& f) {
	File* output_file = new File(schema_, "", Flags::TEMP);
	auto writer = [output_file](TextBuffer* buf) { output_file->out_file_->write(buf->data(), buf->size()); };
	ReorderQueue<TextBuffer*, decltype(writer)> queue(0, writer);

	function<void(int64_t chunk, const Table& t)> callback = [&f, &queue](int64_t chunk, const Table& t) {
		TextBuffer* out = new TextBuffer;
		for (int64_t i = 0; i < t.size(); ++i)
			f(t[i]).write(*out);
		queue.push(chunk, out);
	};

	read(threads, callback);

	return output_file;
}

File* File::sort(int column, int threads) {
	static const int64_t READ_SIZE = 100000;
	file_->rewind();
	Table t(schema_);
	vector<File*> files;
	do {
		t = read(READ_SIZE, threads);
		if (t.size() == 0)
			break;
		files.push_back(new File(schema_, "", Flags::TEMP));
		files.back()->write(t.sorted(column, threads));
		files.back()->out_file_->rewind();
		files.back()->file_->rewind();
	} while (t.size() == READ_SIZE);
	File* out = merge(files.begin(), files.end(), column);
	for (File* f : files)
		delete f;
	return out;
}

}}