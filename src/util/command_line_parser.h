/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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
			throw std::runtime_error("Invalid parameter count for option '--" + id + "'");
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
	struct Add_f
	{
		Add_f(OptionsGroup& parent) :
			parent_(&parent)
		{}
		template<typename T>
		Add_f& operator()(const char* id, char short_id, const char* desc, T& store, T def = T(), const int min_count = 1)
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
	Add_f add()
	{
		return Add_f(*this);
	}
	OptionsGroup(const char* title, const std::vector<unsigned>& commands, bool disabled, CommandLineParserBase* parent) :
		title(title),
		commands(commands),
		disabled(disabled),
		parent_(parent)
	{}
	std::list<std::unique_ptr<OptionBase>> options;
	std::string title;
	const std::vector<unsigned> commands;
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
	OptionsGroup& add_group(const char *title, const std::vector<unsigned>& commands, bool disabled = false);

private:
	void store_option(const std::vector<std::string> &v, unsigned command);

	std::list<OptionsGroup> groups_;
	std::map<std::string, unsigned> command_codes_;
	std::vector<std::pair<std::string, std::string>> commands_;
};