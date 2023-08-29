#pragma once
#include <stdint.h>
#include <cstddef>
#include "def.h"
#include "../text_buffer.h"

namespace Util { namespace Tsv {
	
struct Record {

	Record(const Schema& schema, const char* begin, const char* end):
		schema_(schema),
		buf_(begin),
		end_(end)
	{}

	template<typename T>
	T get(int i) const;
	std::string get(int i) const;

	struct Iterator {
		using iterator_category = std::forward_iterator_tag;
		using difference_type = ptrdiff_t;
		using value_type = const char*;
		using pointer = const char**;
		using reference = const char*&;
		Iterator(Schema::const_iterator it, const char* ptr);
		Iterator& operator++();
		std::string operator*() const;
		bool operator!=(const Iterator& it) const {
			return it_ != it.it_;
		}
		Type type() const {
			return *it_;
		}
		template<typename T>
		T get() const;
	private:
		Schema::const_iterator it_;
		const char* ptr_;
	};

	Iterator begin() const {
		return Iterator(schema_.begin(), buf_);
	}

	Iterator end() const {
		return Iterator(schema_.end(), nullptr);
	}

	int64_t raw_size() const {
		return end_ - buf_;
	}

	void write(TextBuffer& buf) const;

private:

	const Schema& schema_;
	const char* buf_, *end_;

	friend struct Table;

};

}}