#pragma once
#include "tsv.h"
#include "def.h"
#include "record.h"

namespace Util { namespace Tsv {

struct Table {
	
	Table(const Schema& schema);
	template<typename It, typename Tok = TokenIterator<std::string::const_iterator, '\t'>>
	Table(const Schema& schema, It begin, It end);
	Table(Table&& table) noexcept;
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
	void sort(int col, int threads);
	Table sorted(int col, int threads);

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

	void push_string(const std::string& s);
	void push_int32(int32_t x);
	void push_int64(int64_t x);
	void push_int64(const std::string& s);

	Schema schema_;
	std::vector<char> data_;
	std::vector<int64_t> limits_;

	friend struct File;

};

}}