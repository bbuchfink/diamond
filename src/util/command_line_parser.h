/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#ifndef COMMAND_LINE_PARSER_H_
#define COMMAND_LINE_PARSER_H_

#include <vector>
#include <map>
#include <string>
#include <assert.h>
#include <stdlib.h>
#include <stdexcept>
#include <stdint.h>
#include "ptr_vector.h"

template<typename _t>
inline void read_option(_t &dst, const std::vector<std::string> &v)
{ assert(false); }

template<>
inline void read_option<int>(int &dst, const std::vector<std::string> &v)
{ dst = atoi(v[0].c_str()); }

template<>
inline void read_option<unsigned>(unsigned &dst, const std::vector<std::string> &v)
{ dst = atoi(v[0].c_str()); }

template<>
inline void read_option<uint64_t>(uint64_t &dst, const std::vector<std::string> &v)
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
inline void read_option<std::vector<std::string> >(std::vector<std::string> &dst, const std::vector<std::string> &v)
{
	dst = v;
}

template<>
inline void read_option<bool>(bool &dst, const std::vector<std::string> &v)
{ dst = true; }

template<typename _t>
inline bool check_pcount(const std::vector<std::string> &v)
{
	return v.size() == 1;
}

template<>
inline bool check_pcount<bool>(const std::vector<std::string> &v)
{
	return v.size() == 0;
}

template<>
inline bool check_pcount<std::vector<std::string> >(const std::vector<std::string> &v)
{
	return v.size() > 0;
}

struct Option_base
{
	Option_base(const std::string &id, char short_id, const std::string &desc) :
		id(id),
		desc(desc),
		short_id (short_id)		
	{}
	virtual void read(const std::vector<std::string> &v) = 0;
	virtual void set_default() = 0;
	virtual ~Option_base()
	{}
	const std::string id, desc;
	const char short_id;
};

template<typename _t>
struct Option : public Option_base
{
	Option(const char *id, char short_id, const char *desc, _t &store, _t def) :
		Option_base(id, short_id, desc),
		default_(def),
		store_(store)
	{ }
	virtual void read(const std::vector<std::string> &v)
	{
		if (!check_pcount<_t>(v))
			throw std::runtime_error("Invalid parameter count for option '--" + id + "'");
		read_option(store_, v);
	}
	virtual void set_default()
	{
		store_ = default_;
	}
	virtual ~Option()
	{}
private:
	const _t default_;
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
		Add_f& operator()(const char *id, char short_id, const char *desc, _t &store, _t def = _t())
		{
			parent_.options.push_back(new Option<_t>(id, short_id, desc, store, def));
			return *this;
		}
	private:
		Options_group &parent_;
	};
	Add_f add()
	{
		return Add_f(*this);
	}
	Options_group(const char *title):
		title (title)
	{}
	PtrVector<Option_base> options;
	std::string title;
};

struct Command_line_parser
{
	Command_line_parser& add(const Options_group &group);
	Command_line_parser& add_command(const char *s, const char *desc);
	void store(int count, const char **str, unsigned &command);
	void print_help();
private:
	void store_option(const std::vector<std::string> &v);

	std::map<std::string, Option_base*> map_;
	std::map<char, Option_base*> map_short_;
	std::vector<const Options_group*> groups_;
	std::vector<std::pair<std::string,std::string> > commands_;
};

#endif
