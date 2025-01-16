#pragma once
#include <vector>
#include <string>
#include <utility>
#include <ostream>
#include <iomanip>
#include <math.h>

namespace Util {

struct Table {

	Table():
		max_len_(0)
	{}

	Table& operator()(const std::string& s, const std::string& t) {
		data_.emplace_back(s, t);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	Table& operator()(const std::string& s, long long n, const char* unit = "") {
		data_.emplace_back(s, std::to_string(n) + unit);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	Table& operator()(const std::string& s, unsigned long long n, const char* unit = "") {
		data_.emplace_back(s, std::to_string(n) + unit);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	Table& operator()(const std::string& s, long n, const char* unit = "") {
		data_.emplace_back(s, std::to_string(n) + unit);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	Table& operator()(const std::string& s, unsigned long n, const char* unit = "") {
		data_.emplace_back(s, std::to_string(n) + unit);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	Table& operator()(const std::string& s, int n, const char* unit = "") {
		data_.emplace_back(s, std::to_string(n) + unit);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	Table& operator()(const std::string& s, unsigned int n, const char* unit = "") {
		data_.emplace_back(s, std::to_string(n) + unit);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	Table& operator()(const std::string& s, double n, const char* unit = "") {
		data_.emplace_back(s, (n >= 100.0 ? std::to_string((int64_t)round(n)) : std::to_string(n)) + unit);
		max_len_ = std::max(max_len_, s.length());
		return *this;
	}

	friend std::ostream& operator<<(std::ostream& str, const Table& table) {
		for (const auto& l : table.data_)
			str << std::setw(table.max_len_) << l.first << "  " << l.second << std::endl;
		return str;
	}

private:

	std::vector<std::pair<std::string, std::string>> data_;
	size_t max_len_;

};

}
