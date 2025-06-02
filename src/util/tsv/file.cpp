/****
DIAMOND protein aligner
Copyright (C) 2019-2024 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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

#include "file.h"
#include "../enum.h"
#include "../data_structures/reorder_queue.h"
#include "build_helper.h"

using std::vector;
using std::function;
using std::string;
using std::runtime_error;
using std::unique_ptr;
using std::pair;
using Util::String::TokenizerBase;
using Util::String::CharTokenizer;

#ifdef _WIN32
static const char* DEFAULT_LINE_DELIMITER = "\r\n";
#else
static const char* DEFAULT_LINE_DELIMITER = "\n";
#endif

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
	record_id_(0)
{
	//const char* input_file_delim = config.line_delimiter.empty() ? "\n" : config_.line_delimiter.c_str();
	if (flag_all(flags, Flags::WRITE | Flags::OVERWRITE))
		throw std::runtime_error("Invalid File flags");
	if (flag_any(flags_, Flags::RECORD_ID_COLUMN) && schema.front() != Type::INT64)
		throw runtime_error("Schema does not contain record_id column.");
	if (flag_any(flags_, Flags::WRITE | Flags::READ_WRITE)) {
		if (flag_any(flags_, Flags::TEMP)) {
			out_file_.reset(new TempFile);
		}
		else {
			const char* mode = flag_any(flags_, Flags::OVERWRITE) ? "w+b" :
				(flag_any(flags_, Flags::READ_WRITE) ? "r+b" : "wb");
			out_file_.reset(new OutputFile(file_name, Compressor::NONE, mode));
		}
	}
	if (flag_any(flags_, Flags::READ_WRITE)) {
		if (flag_any(flags_, Flags::TEMP)) {
			file_.reset(new TextInputFile(*(TempFile*)out_file_.get()));
		}
		else {
			file_.reset(new TextInputFile(*out_file_));
		}
	}
	else if (!flag_any(flags_, Flags::WRITE))
		file_.reset(new TextInputFile(file_name));
	/*if (config_.line_delimiter.empty())
		config_.line_delimiter = flag_any(flags, Flags::WRITE | Flags::READ_WRITE | Flags::OVERWRITE) ? DEFAULT_LINE_DELIMITER : "\n";*/
}

File::File(const Schema& schema, const std::string& file_name, Flags flags, const Config& config):
	File(schema, file_name.c_str(), flags, config)
{}

File::File(const Schema& schema, unique_ptr<TextInputFile>&& file, Flags flags, const Config& config):
	flags_(get_flags(flags)),
	schema_(schema),
	config_(config),
	record_id_(0)
{
	if (flags != Flags::READ)
		throw runtime_error("File::File");
	file_ = std::move(file);
	//file_->set_separator(config_.line_delimiter.c_str());
}

void File::rewind() {
	if (out_file_)
		out_file_->rewind();
	if (file_)
		file_->rewind();
	record_id_ = 0;
}

bool File::eof() {
	return file_->eof();
}

int64_t File::size() {
	return file_->file_size();
}

File::~File() {
	close();
}

bool File::eof() const {
	return file_->eof();
}

std::string File::file_name() const {
	return file_->file_name;
}

void File::close() {
	if (file_ && file_->temp_file)
		file_->close_and_delete();
	else {
		if (out_file_)
			out_file_->close();
		else if(file_)
			file_->close();
	}
	file_.reset();
	out_file_.reset();
}

void File::read(int threads, std::function<void(int64_t chunk, const Table&)>& callback) {
	function<void(int64_t, const char*, const char*)> f([this, &callback](int64_t chunk, const char* begin, const char* end) {
		Table table(schema_);
		table.append(begin, end, config_.line_tokenizer.get());
		callback(chunk, table);
		});
	read(INT64_MAX, threads, f);
}

Table File::read(int64_t max_size, int threads) {
	BuildHelper helper(std::min(max_size, file_->file_size()));
	auto callback = std::bind(&BuildHelper::add, &helper, std::placeholders::_1);
	ReorderQueue<Table*, decltype(callback)> queue(0, callback);
	function<void(int64_t, const char*, const char*)> f([this, &queue](int64_t chunk, const char* begin, const char* end) {
		Table* table = new Table(this->schema());
		unique_ptr<TokenizerBase> tok(config_.line_tokenizer->clone());
		tok->reset(begin, end);
		pair<bool, string> record;
		while ((record = tok->read_record()).first) {
			table->push_back(record.second.c_str(), record.second.c_str() + record.second.length(), config_.line_tokenizer.get());
		}
		queue.push(chunk, table);
	});
	read(max_size, threads, f);
	return helper.get(schema_);
}

Table File::read(int threads) {
	rewind();
	return read(INT64_MAX, threads);
}

Table File::read_record() {
	Table table(schema_);
	pair<bool, string> r = config_.line_tokenizer->read_record(*file_);
	if(!r.first)
		return table;
	table.push_back(r.second.c_str(), r.second.c_str() + r.second.length(), config_.line_tokenizer.get(), flag_any(flags_, Flags::RECORD_ID_COLUMN) ? record_id_ : -1);
	++record_id_;
	return table;
}

void File::write_record(int i) {
	if (i != (flag_any(flags_, Flags::RECORD_ID_COLUMN) ? 1 : 0))
		throw runtime_error("write_record with insufficient field count.");
	write_buf_ << '\n'; // config_.line_delimiter;
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

void File::write(const TextBuffer* buf) {
	out_file_->write_raw(buf->data(), buf->size());
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
	static const int64_t READ_SIZE = 1000000;
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