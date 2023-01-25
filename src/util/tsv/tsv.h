#pragma once

#include <string>
#include <functional>
#include "../io/text_input_file.h"

namespace Util { namespace Tsv {

std::string fetch_block(TextInputFile& f, std::string& buf);
std::string column(const std::string& line, const size_t i);
std::string columns(const std::string& line, const size_t begin, const size_t end);
size_t column_count(const std::string& line);
std::vector<std::string> extract_column(const std::string& buf, const size_t i);
int64_t count_lines(const std::string& file_name);
template<typename T>
T convert_string(const char* s);
template<typename T>
T convert_string(const std::string& s) {
	return convert_string<T>(s.c_str());
}

template<typename It, char delimiter>
struct TokenIterator {
	TokenIterator(It begin, It end):
		ptr_(begin),
		end_(end)
	{}
	bool good() const {
		return ptr_ < end_;
	}
	std::string operator*() const {
		return std::string(ptr_, std::find(ptr_, end_, delimiter));
	}
	TokenIterator& operator++() {
		ptr_ = std::find(ptr_, end_, delimiter);
#ifndef NDEBUG
		if (ptr_ < end_)
#endif
			++ptr_;
		return *this;
	}
	It ptr() const {
		return ptr_;
	}
private:
	It ptr_, end_;
};

using LineIterator = TokenIterator<const char*, '\n'>;

template<typename It>
struct FastaIterator {
	FastaIterator(It begin, It end):
		ptr_(begin),
		end_(end)
	{}
private:
	It ptr_, end_;
};

}}