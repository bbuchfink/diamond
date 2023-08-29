#pragma once
#include "tsv.h"
#include "def.h"
#include "record.h"

namespace Util { namespace Tsv {

struct File;
struct Table;

using MapFunc = std::function<Table(const Record&)>;

struct Table {
	
	Table(const Schema& schema);
	template<typename It, typename Tok = TokenIterator<std::string::const_iterator, '\t'>>
	Table(const Schema& schema, It begin, It end);
	Table(Table&& table) noexcept;
	Table(const Schema& schema, std::vector<char>&& data, std::vector<int64_t>&& limits);
	Table& operator=(Table&& table) noexcept;

	const Schema& schema() const {
		return schema_;
	}

	int64_t size() const {
		return limits_.size() - 1;
	}

	bool empty() const {
		return size() == 0;
	}

	Record operator[](int64_t i) const {
		return Record(schema_, data_.data() + limits_[i], data_.data() + limits_[i + 1]);
	}

	Record front() const {
		return Record(schema_, data_.data(), data_.data() + limits_[1]);
	}

	template<typename It, typename Tok = TokenIterator<It, '\t'>>
	void push_back(It begin, It end, RecordId record_id = -1);
	void push_back(const Record& record);
	void append(const Table& table);
	template<typename It, typename Tok = TokenIterator<std::string::const_iterator, '\t'>>
	void append(It begin, It end);
	void write(TextBuffer& buf) const;

	template<typename... Targs>
	void write_record(Targs... FArgs) {
		limits_.push_back(limits_.back());
		write_record((int)schema_.size(), FArgs...);
	};

	void sort(int col, int threads);
	void map(int threads, MapFunc& f, File& out) const;
	Table sorted(int col, int threads);
	int64_t alloc_size() const;

	template<typename It>
	Table shuffle(It begin, It end) {
		Table t(schema_);
		t.data_.reserve(data_.size());
		t.limits_.reserve(limits_.size());
		for (It i = begin; i != end; ++i)
			t.push_back(operator[](*i));
		return t;
	}
	
private:

	void write_record(int i);

	template<typename T, typename... Targs>
	void write_record(int i, T value, Targs... FArgs) {
		if (i == 0)
			throw std::runtime_error("write_record with too many fields.");
		push(value);
		write_record(i - 1, FArgs...);
	}

	void push(const std::string& s);
	void push(int32_t x);
	void push(int64_t x);
	
	Schema schema_;
	std::vector<char> data_;
	std::vector<int64_t> limits_;

	friend struct File;
	friend struct BuildHelper;

};

}}