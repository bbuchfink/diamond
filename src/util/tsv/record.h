/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once
#include <stdint.h>
#include <stddef.h>
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