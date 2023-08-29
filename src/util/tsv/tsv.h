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
	bool good() const {
		return ptr_ < end_;
	}
	std::string operator*() const {
		if (*ptr_ == '>') {
			It i = std::find(ptr_, end_, '\n');
			return std::string(ptr_ + 1, i);
		}
		else {
			std::string r;
			r.reserve(end_ - ptr_);
			It p = ptr_;
			while (p < end_) {
				It i = std::find(p, end_, '\n');
				r.append(p, i);
				p = i + 1;
			}
			return r;
		}
	}
	FastaIterator& operator++() {
		if (*ptr_ == '>') {
			ptr_ = std::find(ptr_, end_, '\n') + 1;
		}
		else {
			throw std::runtime_error("Seeking FASTA iterator past end.");
		}
		return *this;
	}
private:
	It ptr_, end_;
};

struct MalformedFastqRecord : public std::exception {
};

template<typename It>
struct FastqIterator {
	FastqIterator(It begin, It end) :
		ptr_(begin),
		end_(end)
	{}
	bool good() const {
		return ptr_ < end_;
	}
	std::string operator*() const {
		if (*ptr_ == '@') {
			It i = std::find(ptr_, end_, '\n');
			return std::string(ptr_ + 1, i);
		}
		else if (*ptr_ == '+') {
			It i = std::find(ptr_, end_, '\n');
			if (i == end_)
				return std::string();
			return std::string(i + 1, end_);
		}
		else {
			return std::string(ptr_, std::find(ptr_, end_, '\n'));
		}
	}
	FastqIterator& operator++() {
		if (*ptr_ == '@') {
			ptr_ = std::find(ptr_, end_, '\n') + 1;
		}
		else if (*ptr_ == '+') {
			throw std::runtime_error("Seeking FASTQ iterator past end.");
		}
		else {
			ptr_ = std::find(ptr_, end_, '\n') + 1;
			if (ptr_ < end_ && *ptr_ != '+')
				throw MalformedFastqRecord();
		}
		return *this;
	}
private:
	It ptr_, end_;
};

}}