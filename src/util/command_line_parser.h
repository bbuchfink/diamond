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

#pragma once
#include <memory>
#include <vector>
#include <list>
#include <map>
#include <string>
#include <assert.h>
#include <stdlib.h>
#include <stdexcept>
#include <stdint.h>
#include "options/option.h"

template<typename T>
inline void read_option(T &dst, const std::vector<std::string> &v)
{
	throw std::runtime_error("Invalid option type.");
}

template<>
inline void read_option<int>(int &dst, const std::vector<std::string> &v)
{ dst = atoi(v[0].c_str()); }

template<>
inline void read_option<unsigned>(unsigned &dst, const std::vector<std::string> &v)
{ dst = atoi(v[0].c_str()); }

template<>
inline void read_option<long>(long& dst, const std::vector<std::string>& v)
{
	dst = (long)atoll(v[0].c_str());
}

template<>
inline void read_option<unsigned long>(unsigned long& dst, const std::vector<std::string>& v)
{
	dst = (unsigned long)atoll(v[0].c_str());
}

template<>
inline void read_option<long long>(long long& dst, const std::vector<std::string>& v)
{
	dst = atoll(v[0].c_str());
}

template<>
inline void read_option<Option<int64_t>>(Option<int64_t>& dst, const std::vector<std::string>& v)
{
	dst = atoll(v[0].c_str());
}

template<>
inline void read_option<unsigned long long>(unsigned long long& dst, const std::vector<std::string>& v)
{
	dst = atoll(v[0].c_str());
}

template<>
inline void read_option<double>(double &dst, const std::vector<std::string> &v)
{
	dst = atof(v[0].c_str());
}


template<>
inline void read_option<Option<double>>(Option<double> &dst, const std::vector<std::string> &v)
{
	dst = atof(v[0].c_str());
}

template<>
inline void read_option<std::string>(std::string &dst, const std::vector<std::string> &v)
{ dst = v[0]; }

template<>
inline void read_option<Option<std::string>>(Option<std::string> &dst, const std::vector<std::string> &v)
{
	dst = v[0];
}

template<>
inline void read_option<std::vector<std::string>>(std::vector<std::string> &dst, const std::vector<std::string> &v)
{
	dst = v;
}

template<>
inline void read_option<Option<std::vector<std::string>>>(Option<std::vector<std::string>>& dst, const std::vector<std::string>& v)
{
	dst = v;
}

template<>
inline void read_option<bool>(bool &dst, const std::vector<std::string> &v)
{ dst = true; }

template<typename T>
inline bool check_pcount(const std::vector<std::string> &v, const int min_count)
{
	return v.size() == 1;
}

template<>
inline bool check_pcount<bool>(const std::vector<std::string> &v, const int min_count)
{
	return v.size() == 0;
}

template<>
inline bool check_pcount<std::vector<std::string>>(const std::vector<std::string> &v, const int min_count)
{
	return v.size() >= (size_t)min_count;
}

template<>
inline bool check_pcount<Option<std::vector<std::string>>>(const std::vector<std::string>& v, const int min_count)
{
	return v.size() >= (size_t)min_count;
}

template<typename T>
inline bool check_present(const T& v) {
	return false;
}

template<>
inline bool check_present<std::string>(const std::string& v) {
	return !v.empty();
}

template<typename T>
inline void set_base_ptr(T& v, const OptionBase* ptr) {
}

template<typename T>
inline void set_base_ptr(Option<T>& v, const OptionBase* ptr) {
	v.base_ = ptr;
}

template<typename T>
inline void set_option_default(T& option, const T& value) {
	option = value;
}

template<typename T>
inline void set_option_default(Option<T>& option, const Option<T>& value) {
}

template<typename T>
struct OptionDesc : public OptionBase
{
	OptionDesc(const char *id, char short_id, const char *desc, bool disabled, T &store, T def, const int min_count, const OptionsGroup* group) :
		OptionBase(id, short_id, desc, disabled, group),
		default_(def),
		min_count(min_count),
		store_(&store)
	{ }
	virtual void read(const std::vector<std::string> &v)
	{
		if (!check_pcount<T>(v, min_count))
			throw std::runtime_error("Invalid parameter count for option '" + (short_id ? std::string("-") + short_id + '/' : std::string()) + "--" + id + "'");
		read_option(*store_, v);
	}
	virtual bool present() {
		return check_present<T>(*store_);			
	}
	virtual void set_default()
	{
		set_option_default(*store_, default_);
	}
	virtual ~OptionDesc()
	{}
private:
	const T default_;
	const int min_count;
	T* store_;
};

struct CommandLineParserBase
{
	std::map<std::string, OptionBase*> map_;
	std::map<char, OptionBase*> map_short_;
};

struct OptionsGroup
{
	struct AddFunc
	{
		AddFunc(OptionsGroup& parent) :
			parent_(&parent)
		{}
		template<typename T>
		AddFunc& operator()(const char* id, char short_id, const char* desc, T& store, T def = T(), const int min_count = 1)
		{
			auto o = new OptionDesc<T>(id, short_id, desc, parent_->disabled, store, def, min_count, parent_);
			parent_->options.emplace_back(o);
			parent_->parent_->map_[id] = o;
			parent_->parent_->map_short_[short_id] = o;
			set_base_ptr(store, o);
			return *this;
		}
	private:
		OptionsGroup* parent_;
	};
	AddFunc add()
	{
		return AddFunc(*this);
	}
	OptionsGroup(const char* title, const std::vector<int>& commands, bool disabled, CommandLineParserBase* parent) :
		title(title),
		commands(commands),
		disabled(disabled),
		parent_(parent)
	{}
	std::list<std::unique_ptr<OptionBase>> options;
	std::string title;
	const std::vector<int> commands;
	bool disabled;
private:
	CommandLineParserBase* parent_;
};

struct CommandLineParser : CommandLineParserBase
{
	CommandLineParser& add_command(const char *s, const char *desc, unsigned code);
	void store(int count, const char **str, unsigned &command);
	void require(const char* option);
	void print_help();
    void print_documentation(int command);
	OptionsGroup& add_group(const char *title, const std::vector<int>& commands, bool disabled = false);

private:
	void store_option(const std::vector<std::string> &v, unsigned command);

	std::list<OptionsGroup> groups_;
	std::map<std::string, unsigned> command_codes_;
	std::vector<std::pair<std::string, std::string>> commands_;
};