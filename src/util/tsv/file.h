/****
DIAMOND protein aligner
Copyright (C) 2019-2024 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink

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

#pragma once
#include <functional>
#include <memory>
#include "../io/text_input_file.h"
#include "table.h"
#include "../enum.h"

namespace Util { namespace Tsv {

enum class Flags {
	READ = 0, READ_WRITE = 1, WRITE = 1 << 1, OVERWRITE = 1 << 2, RECORD_ID_COLUMN = 1 << 3, TEMP = 1 << 4
};

DEFINE_ENUM_FLAG_OPERATORS(Flags);

struct FileColumn {
	int file, column;
};

struct Config {
	Config(const char* line_delimiter = "", Util::String::TokenizerBase* line_tokenizer = new Util::String::CharTokenizer('\t')) :
		line_delimiter(line_delimiter),
		line_tokenizer(line_tokenizer)
	{}
	std::string line_delimiter;
	std::shared_ptr<Util::String::TokenizerBase> file_tokenizer, line_tokenizer;
};

struct File {

	File(const Schema& schema, const char* file_name, Flags flags = Flags(), const Config& config = Config());
	File(const Schema& schema, const std::string& file_name, Flags flags = Flags(), const Config& config = Config());
	File(const Schema& schema, std::unique_ptr<TextInputFile>&& file, Flags flags = Flags(), const Config& config = Config());
	void rewind();
	bool eof();
	int64_t size();
	~File();

	void read(int64_t max_size, int threads, std::function<void(int64_t chunk, const char*, const char*)>& callback);
	void read(int threads, std::function<void(int64_t chunk, const Table&)>& callback);
	Table read(int64_t max_size, int threads);
	Table read(int threads);
	Table read_record();

	template<typename... Targs, typename Out> void read(Out out);

	template<typename F, typename Out>
	void read(int threads, F& f, Out& out) {
		/*using T = typename std::result_of<F(const Record&)>::type;
		auto f2 = [](std::vector<T>* v) {
			for (std::vector<T>::const_iterator i = v->cbegin(); i != v->cend(); ++i)
				*out++ = *i;
		}
		ReorderQueue<std::vector<T>*, decltype(f2)> queue;
		auto f1 = [&queue](int64_t chunk const Table& table) {
			std::vector<T>* v = new vector<T>();
			const int64_t n = table.size();
			v->reserve(n);
			for (int64_t i = 0; i < n; ++i)
				v->push_back(f(table[i]));
			queue.push(chunk, v);
		};*/
	}

	Schema schema() const {
		return schema_;
	}

	void write(const Record& record);
	void write(const Table& table);
	void write(const TextBuffer* buf);

	template<typename... Targs>
	void write_record(Targs... FArgs) {
		write_record((int)schema_.size(), FArgs...);
	}

	/*template<typename... Targs>
	void write_record(const std::tuple<Targs...>& tuple) {
		TupleDispatch<Targs...>()(*this, tuple, typename gens<sizeof...(Targs)>::type());
	}*/

	template<typename It, typename F>
	void write(It begin, It end, F& f) {
		for (It i = begin; i != end; ++i)
			write_record(f(*i));
	}

	File* map(int threads, std::function<Table(const Record&)>& f);

	File* sort(int column, int threads);
	void close();
	bool eof() const;
	std::string file_name() const;
	
private:

	void write_record(int i);

	template<typename T, typename... Targs>
	void write_record(int i, T value, Targs... FArgs) {
		if (i == 0 || (i == 1 && flag_any(flags_, Flags::RECORD_ID_COLUMN)))
			throw std::runtime_error("write_record with too many fields.");
		write_buf_ << value;
		if (i > 1)
			write_buf_ << '\t';
		write_record(i - 1, FArgs...);
	}

	/*template<typename... Targs>
	struct TupleDispatch {
		template<int ...S>
		void operator()(File& f, const std::tuple<Targs...>& tuple, seq<S...>) const {
			f.write_record((int)f.schema_.size(), std::get<S>(tuple) ...);
		}
	};*/

	const Flags flags_;
	const Schema schema_;
	Config config_;
	std::unique_ptr<OutputFile> out_file_;
	std::unique_ptr<TextInputFile> file_;
	TextBuffer write_buf_;
	RecordId record_id_;

	friend File* merge(std::vector<File*>::iterator begin, std::vector<File*>::iterator end, int column);
	friend void join(File& file1, File& file2, int column1, int column2, const std::vector<FileColumn>& output_fields, File& out);
	friend File* join(File& file1, File& file2, int column1, int column2, const std::vector<FileColumn>& output_fields);

};

}}