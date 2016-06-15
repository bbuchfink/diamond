/****
Copyright (c) 2016, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

using std::vector;
using std::string;
using std::pair;

template<typename _t>
inline void read_option(_t &dst, const vector<string> &v)
{ assert(false); }

template<>
inline void read_option<int>(int &dst, const vector<string> &v)
{ dst = atoi(v[0].c_str()); }

template<>
inline void read_option<unsigned>(unsigned &dst, const vector<string> &v)
{ dst = atoi(v[0].c_str()); }

template<>
inline void read_option<uint64_t>(uint64_t &dst, const vector<string> &v)
{
	dst = atoi(v[0].c_str());
}

template<>
inline void read_option<double>(double &dst, const vector<string> &v)
{
	dst = atof(v[0].c_str());
}

template<>
inline void read_option<string>(string &dst, const vector<string> &v)
{ dst = v[0]; }

template<>
inline void read_option<bool>(bool &dst, const vector<string> &v)
{ dst = true; }

template<typename _t>
inline bool check_pcount(const vector<string> &v)
{
	return v.size() == 1;
}

template<>
inline bool check_pcount<bool>(const vector<string> &v)
{
	return v.size() == 0;
}

struct Option_base
{
	Option_base(const string &id, char short_id, const string &desc) :
		id(id),
		desc(desc),
		short_id (short_id)		
	{}
	virtual void read(const vector<string> &v) = 0;
	virtual void set_default() = 0;
	virtual ~Option_base()
	{}
	const string id, desc;
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
	virtual void read(const vector<string> &v)
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
	Ptr_vector<Option_base> options;
	string title;
};

struct Command_line_parser
{
	Command_line_parser& add(const Options_group &group);
	Command_line_parser& add_command(const char *s, const char *desc);
	void store(int count, const char **str, unsigned &command);
	void print_help();
private:
	void store_option(const vector<string> &v);

	std::map<string, Option_base*> map_;
	std::map<char, Option_base*> map_short_;
	vector<const Options_group*> groups_;
	vector<pair<string,string> > commands_;
};

#endif