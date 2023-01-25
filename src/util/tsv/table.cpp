#define _REENTRANT
#include "../../lib/ips4o/ips4o.hpp"
#include "table.h"
#include "../sequence/sequence.h"
#include "../algo/sort_helper.h"
#include "../algo/transform_iterator.h"

using std::numeric_limits;
using std::runtime_error;
using std::string;
using std::vector;
using std::pair;
using std::less;
using std::move;

namespace Util { namespace Tsv {

Table::Table(const Schema& schema) :
	schema_(schema),
	limits_({ 0 })
{
}

Table::Table(Table&& table) noexcept:
	schema_(table.schema_),
	data_(move(table.data_)),
	limits_(move(table.limits_))
{}

template<typename It, typename Tok>
Table::Table(const Schema& schema, It begin, It end):
	Table(schema)
{
	append(begin, end);
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

template<typename It, typename Tok>
void Table::append(It begin, It end) {
	LineIterator it(begin, end);
	string line;
	while (it.good()) {
		line = *it;
		push_back(line.cbegin(), line.cend(), -1);
		++it;
	}
}

template void Table::append<const char*, TokenIterator<string::const_iterator, '\t'>>(const char*, const char*);
//template Table::Table<const char*, TokenIterator<string::const_iterator, '\t'>>(const Schema&, const char*, const char*);

Table& Table::operator=(Table&& table) noexcept {
	schema_ = table.schema_;
	limits_ = move(table.limits_);
	data_ = move(table.data_);
	return *this;
}

void Table::push_back(const Record& record) {
	limits_.push_back(limits_.back() + record.raw_size());
	data_.insert(data_.end(), record.buf_, record.end_);
}

template<typename It, typename Tok>
void Table::push_back(It begin, It end, RecordId record_id) {
	Tok tok(begin, end);
	Schema::const_iterator i = schema_.cbegin();
	limits_.push_back(limits_.back());
	if (record_id >= 0) {
		++i;
		push_int64(record_id);
	}
	while (tok.good() && i < schema_.cend()) {
		switch (*i) {
		case Type::STRING:
			push_string(*tok);
			break;
		case Type::INT64:
			push_int64(*tok);
			break;
		default:
			throw runtime_error("Invalid type in schema");
		}
		++tok;
		++i;
	}
	if (i < schema_.cend())
		throw runtime_error("Missing fields in input line");
}

template void Table::push_back<string::const_iterator, TokenIterator<string::const_iterator, '\t'>>(string::const_iterator, string::const_iterator, RecordId);
template void Table::push_back<const char*, LineIterator>(const char*, const char*, RecordId);

void Table::push_string(const std::string& s) {
	push_int32((int32_t)s.length());
	data_.insert(data_.end(), s.cbegin(), s.cend());
	limits_.back() += s.length();
}

void Table::push_int32(int32_t x) {
	data_.insert(data_.end(), reinterpret_cast<const char*>(&x), reinterpret_cast<const char*>(&x) + 4);
	limits_.back() += 4;
}

void Table::push_int64(int64_t x) {
	data_.insert(data_.end(), reinterpret_cast<const char*>(&x), reinterpret_cast<const char*>(&x) + 8);
	limits_.back() += 8;
}

void Table::push_int64(const std::string& s) {
	push_int64(convert_string<int64_t>(s));
}

void Table::write(TextBuffer& buf) const {
	for (int64_t i = 0; i < size(); ++i)
		(*this)[i].write(buf);
}

Table Table::sorted(int col, int threads) {
	if ((size_t)col >= schema_.size() || schema_[col] != Type::INT64)
		throw runtime_error("Invalid sort");
	vector<pair<int64_t, int64_t>> v;
	v.reserve(size());
	for (int64_t i = 0; i < size(); ++i)
		v.emplace_back(this->operator[](i).get<int64_t>(col), i);
	ips4o::parallel::sort(v.begin(), v.end(), less<>(), threads);
	return shuffle(transform(v.begin(), Second<int64_t, int64_t>()), transform(v.end(), Second<int64_t, int64_t>()));
}

void Table::sort(int col, int threads) {
	*this = sorted(col, threads);
}

}}