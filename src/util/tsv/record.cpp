#include <string>
#include "record.h"

using std::advance;
using std::string;
using std::to_string;

namespace Util { namespace Tsv {

static int64_t interpret(const char* ptr, int64_t) {
	return *reinterpret_cast<const int64_t*>(ptr);
}

static int32_t interpret(const char* ptr, int32_t) {
	return *reinterpret_cast<const int32_t*>(ptr);
}

static string interpret(const char* ptr, const string&) {
	return string(ptr + 4, ptr + 4 + *reinterpret_cast<const int32_t*>(ptr));
}

Record::Iterator::Iterator(Schema::const_iterator it, const char* ptr):
	it_(it),
	ptr_(ptr)
{}

Record::Iterator& Record::Iterator::operator++() {
	switch (*it_) {
	case Type::STRING:
		ptr_ += 4 + *(int32_t*)ptr_;
		break;
	case Type::INT64:
		ptr_ += 8;
		break;
	}
	++it_;
	return *this;
}

string Record::Iterator::operator*() const {
	switch (*it_) {
	case Type::STRING:
		return interpret(ptr_, string());
	case Type::INT64:
		return to_string(interpret(ptr_, int64_t()));
	default:
		throw InvalidType();
	}
}

template<typename T>
T Record::Iterator::get() const {
	return interpret(ptr_, T());
}

template int64_t Record::Iterator::get<int64_t>() const;
template string Record::Iterator::get<string>() const;

template<typename T>
T Record::get(int i) const {
	Iterator it = begin();
	advance(it, i);
	return it.get<T>();
}

template int64_t Record::get<int64_t>(int i) const;
template string Record::get<string>(int i) const;

string Record::get(int i) const {
	Iterator it = begin();
	advance(it, i);
	return *it;
}

void Record::write(TextBuffer& buf) const {
	int n = 0;
	for (Record::Iterator i = begin(); i != end(); ++i) {
		if (n++ > 0)
			buf << '\t';
		switch (i.type()) {
		case Type::INT64:
			buf << i.get<int64_t>();
			break;
		case Type::STRING:
			buf << i.get<string>();
			break;
		}
	}
	buf << '\n';
}

}}