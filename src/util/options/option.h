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

#include <string>
#include <vector>
#include <stdexcept>
#pragma once

struct OptionsGroup;

struct OptionBase
{
	OptionBase(const std::string& id, char short_id, const std::string& desc, bool disabled, const OptionsGroup* group) :
		id(id),
		desc(desc),
		short_id(short_id),
		disabled(disabled),
		group(group)
	{}
	virtual void read(const std::vector<std::string>& v) = 0;
	virtual bool present() = 0;
	virtual void set_default() = 0;
	virtual ~OptionBase()
	{}
	const std::string id, desc;
	const char short_id;
	bool disabled;
	const OptionsGroup* group;
};

template<typename T>
struct Option : public T {
	Option():
		present_(false)
	{}
	bool present() const {
		return present_;
	}
	bool blank() const {
		return !present_;
	}
	Option& operator=(const T& value) {
		*static_cast<T*>(this) = value;
		present_ = true;
		return *this;
	}
	void require() {
		if(!present_)
			throw std::runtime_error("Missing parameter: --" + base_->id + "/-" + base_->short_id);
	}
	T get(const T& default_value) const {
		return present_ ? *this : default_value;
	}
	T get_present() const {
		if (!present_)
			throw std::runtime_error("Option::get_present");
		return *this;
	}
	void unset() {
		present_ = false;
	}
private:
	bool present_;
	const OptionBase* base_;
	template<typename T2> friend void set_base_ptr(Option<T2>&, const OptionBase*);
};

template<>
struct Option<double> {
	Option() :
		present_(false)
	{}
	bool present() const {
		return present_;
	}
	bool blank() const {
		return !present_;
	}
	Option& operator=(const double value) {
		value_ = value;
		present_ = true;
		return *this;
	}
	void set_if_blank(const double value) {
		if (!present_)
			this->operator=(value);
	}
	void unset() {
		present_ = false;
	}
	operator double() const {
		if (!present_)
			throw std::runtime_error("Option::present");
		return value_;
	}
	double get(const double default_value) const {
		return present_ ? value_ : default_value;
	}
	double get_present() const {
		if (!present_)
			throw std::runtime_error("Option::get_present");
		return value_;
	}
private:
	double value_;
	bool present_;
	const OptionBase* base_;
	template<typename T2> friend void set_base_ptr(Option<T2>&, const OptionBase*);
};

template<>
struct Option<int64_t> {
	Option() :
		present_(false)
	{}
	bool present() const {
		return present_;
	}
	bool blank() const {
		return !present_;
	}
	Option& operator=(const int64_t value) {
		value_ = value;
		present_ = true;
		return *this;
	}
	void set_if_blank(const int64_t value) {
		if (!present_)
			this->operator=(value);
	}
	void unset() {
		present_ = false;
	}
	operator int64_t() const {
		if (!present_)
			throw std::runtime_error("Option::present");
		return value_;
	}
	int64_t get(const int64_t default_value) const {
		return present_ ? value_ : default_value;
	}
private:
	int64_t value_;
	bool present_;
	const OptionBase* base_;
	template<typename T2> friend void set_base_ptr(Option<T2>&, const OptionBase*);
};