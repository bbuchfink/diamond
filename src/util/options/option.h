/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

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