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
#include <vector>
#include <map>
#include <string>
#include <assert.h>
#include <stdlib.h>
#include <stdexcept>
#include <stdint.h>
#include "ptr_vector.h"
#include "options/option.h"

template<typename _t>
inline void read_option(_t &dst, const std::vector<std::string> &v)
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
inline void read_option<std::string>(std::string &dst, const std::vector<std::string> &v)
{ dst = v[0]; }

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

template<typename _t>
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


struct Option_base
{
	Option_base(const std::string &id, char short_id, const std::string &desc, bool disabled) :
		id(id),
		desc(desc),
		short_id (short_id),
		disabled(disabled)
	{}
	virtual void read(const std::vector<std::string> &v) = 0;
	virtual void set_default() = 0;
	virtual ~Option_base()
	{}
	const std::string id, desc;
	const char short_id;
	bool disabled;
};

template<typename _t>
struct OptionDesc : public Option_base
{
	OptionDesc(const char *id, char short_id, const char *desc, bool disabled, _t &store, _t def, const int min_count) :
		Option_base(id, short_id, desc, disabled),
		default_(def),
		min_count(min_count),
		store_(store)
	{ }
	virtual void read(const std::vector<std::string> &v)
	{
		if (!check_pcount<_t>(v, min_count))
			throw std::runtime_error("Invalid parameter count for option '--" + id + "'");
		read_option(store_, v);
	}
	virtual void set_default()
	{
		store_ = default_;
	}
	virtual ~OptionDesc()
	{}
private:
	const _t default_;
	const int min_count;
	_t &store_;
};

struct Options_group
{
	struct Add_f
	{
		Add_f(Options_group &parent):
			parent_(parent)
		{}
		template<typename _t>
		Add_f& operator()(const char *id, char short_id, const char *desc, _t &store, _t def = _t(), const int min_count = 1)
		{
			parent_.options.push_back(new OptionDesc<_t>(id, short_id, desc, parent_.disabled, store, def, min_count));
			return *this;
		}
	private:
		Options_group &parent_;
	};
	Add_f add()
	{
		return Add_f(*this);
	}
	Options_group(const char *title, bool disabled = false):
		title (title),
		disabled(disabled)
	{}
	PtrVector<Option_base> options;
	std::string title;
	bool disabled;
};

struct Command_line_parser
{
	Command_line_parser& add(const Options_group &group);
	Command_line_parser& add_command(const char *s, const char *desc, unsigned code);
	void store(int count, const char **str, unsigned &command);
	void print_help();
private:
	void store_option(const std::vector<std::string> &v);

	std::map<std::string, Option_base*> map_;
	std::map<char, Option_base*> map_short_;
	std::vector<const Options_group*> groups_;
	std::map<std::string, unsigned> command_codes_;
	std::vector<std::pair<std::string, std::string>> commands_;
};