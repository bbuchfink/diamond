#define _REENTRANT
#include "ips4o/ips4o.hpp"
#include "table.h"
#include "../algo/sort_helper.h"
#include "../algo/transform_iterator.h"
#include "file.h"
#include "../data_structures/reorder_queue.h"
#include "../util.h"

using std::numeric_limits;
using std::runtime_error;
using std::string;
using std::vector;
using std::pair;
using std::less;
using std::bind;
using std::atomic;
using std::thread;
using std::unique_ptr;
using Util::String::TokenizerBase;

namespace Util { namespace Tsv {

Table::Table(const Schema& schema) :
	schema_(schema),
	limits_({ 0 })
{
}

Table::Table(Table&& table) noexcept:
	schema_(table.schema_),
	data_(std::move(table.data_)),
	limits_(std::move(table.limits_))
{}

Table::Table(const Schema& schema, std::vector<char>&& data, std::vector<int64_t>&& limits):
	schema_(schema),
	data_(std::move(data)),
	limits_(std::move(limits))
{
}

template<typename It, typename Tok>
Table::Table(const Schema& schema, It begin, It end):
	Table(schema)
{
	append(begin, end);
}

int64_t Table::alloc_size() const {
	return data_.size() + limits_.size() * sizeof(decltype(limits_)::value_type);
}

void Table::append(const Table& table) {
	if (schema_ != table.schema_)
		throw SchemaMismatch();
	data_.insert(data_.end(), table.data_.begin(), table.data_.end());
	const int64_t d = limits_.back();
	limits_.reserve(limits_.size() + table.limits_.size() - 1);
	for (auto it = table.limits_.cbegin() + 1; it != table.limits_.cend(); ++it)
		limits_.push_back(*it + d);
}

void Table::append(const char* begin, const char* end, const TokenizerBase* tok) {
	unique_ptr<TokenizerBase> it(tok->clone());
	it->reset(begin, end);
	string line;
	while (it->good()) {
		line = **it;
		push_back(line.c_str(), line.c_str() + line.length(), tok, -1);
		++(*it);
	}
}

Table& Table::operator=(Table&& table) noexcept {
	schema_ = table.schema_;
	limits_ = std::move(table.limits_);
	data_ = std::move(table.data_);
	return *this;
}

void Table::push_back(const Record& record) {
	limits_.push_back(limits_.back() + record.raw_size());
	data_.insert(data_.end(), record.buf_, record.end_);
}

void Table::push_back(const char* begin, const char* end, const TokenizerBase* tok, RecordId record_id) {
	unique_ptr<TokenizerBase> it(tok->clone());
	it->reset(begin, end);
	Schema::const_iterator i = schema_.cbegin();
	limits_.push_back(limits_.back());
	if (record_id >= 0) {
		++i;
		push(record_id);
	}
	while (it->good() && i < schema_.cend()) {
		switch (*i) {
		case Type::STRING:
			push(**it);
			break;
		case Type::INT64:
			push(Util::String::convert_string<int64_t>(**it));
			break;
		default:
			throw runtime_error("Invalid type in schema");
		}
		++(*it);
		++i;
	}
	if (i < schema_.cend())
		throw runtime_error("Missing fields in input line");
}

void Table::push(const std::string& s) {
	push((int32_t)s.length());
	data_.insert(data_.end(), s.cbegin(), s.cend());
	limits_.back() += s.length();
}

void Table::push(int32_t x) {
	data_.insert(data_.end(), reinterpret_cast<const char*>(&x), reinterpret_cast<const char*>(&x) + 4);
	limits_.back() += 4;
}

void Table::push(int64_t x) {
	data_.insert(data_.end(), reinterpret_cast<const char*>(&x), reinterpret_cast<const char*>(&x) + 8);
	limits_.back() += 8;
}

void Table::write(TextBuffer& buf) const {
	for (int64_t i = 0; i < size(); ++i)
		(*this)[i].write(buf);
}

void Table::write_record(int i) {
	if (i != 0)
		throw std::runtime_error("Mismatching field count for Table::write_record");
}

Table Table::sorted(int col, int threads) {
	if ((size_t)col >= schema_.size() || schema_[col] != Type::INT64)
		throw runtime_error("Invalid sort");
	vector<pair<int64_t, int64_t>> v;
	v.reserve(size());
	for (int64_t i = 0; i < size(); ++i)
		v.emplace_back(this->operator[](i).get<int64_t>(col), i);
	ips4o::parallel::sort(v.begin(), v.end(), less<pair<int64_t, int64_t>>(), threads);
	return shuffle(transform(v.begin(), Second<int64_t, int64_t>()), transform(v.end(), Second<int64_t, int64_t>()));
}

void Table::sort(int col, int threads) {
	*this = sorted(col, threads);
}

void Table::map(int threads, std::function<Table(const Record&)>& f, File& out) const {
	static const RecordId BATCH_SIZE = 1024;
	auto callback = bind(static_cast<void (File::*)(const TextBuffer*)>(&File::write), &out, std::placeholders::_1);
	ReorderQueue<TextBuffer*, decltype(callback)> queue(0, callback);
	atomic<int64_t> next(0);
	auto worker = [&]() {
		for (;;) {
			const int64_t p = next.fetch_add(1, std::memory_order_relaxed);
			if (p * BATCH_SIZE >= size())
				break;
			TextBuffer* buf = new TextBuffer;
			for (RecordId i = p * BATCH_SIZE; i < std::min((p + 1) * BATCH_SIZE, size()); ++i)
				f((*this)[i]).write(*buf);
			queue.push(p, buf);
		}
	};
	vector<thread> t;
	for (int i = 0; i < std::min((int64_t)threads, div_up(size(), BATCH_SIZE)); ++i)
		t.emplace_back(worker);
	for (auto& i : t)
		i.join();
}

}}