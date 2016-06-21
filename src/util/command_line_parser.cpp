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

#include <iostream>
#include <algorithm>
#include "command_line_parser.h"

using std::cout;
using std::endl;

Command_line_parser & Command_line_parser::add(const Options_group & group)
{
	for (Ptr_vector<Option_base>::const_iterator i = group.options.begin(); i != group.options.end(); ++i) {
		map_[(*i)->id] = *i;
		map_short_[(*i)->short_id] = *i;
	}
	groups_.push_back(&group);
	return *this;
}

Command_line_parser & Command_line_parser::add_command(const char * s, const char *desc)
{
	commands_.push_back(pair<string,string>(s,desc));
	return *this;
}

void Command_line_parser::store_option(const vector<string> &v)
{
	if (v.size() == 0)
		return;
	Option_base* o = 0;
	string id;
	vector<string> v2;

	if (v[0].length() <= 1)
		throw std::runtime_error("Invalid option syntax.");
	else if (v[0].substr(0, 2) == "--") {
		id = v[0].substr(2);
		std::map<string, Option_base*>::const_iterator i = map_.find(id);
		if (i != map_.end())
			o = i->second;
	} else if (v[0][0] == '-') {
		id = v[0][1];
		std::map<char, Option_base*>::const_iterator i = map_short_.find(v[0][1]);
		if (i != map_short_.end())
			o = i->second;
		if (v[0].length() > 2)
			v2.push_back(v[0].substr(2));
	} else
		throw std::runtime_error("Command line options must begin with - or --.");
	
	if(o == 0)
		throw std::runtime_error("Invalid option: " + id);
	else {
		v2.insert(v2.end(), v.begin() + 1, v.end());
		o->read(v2);
	}
}

void Command_line_parser::store(int count, const char ** str, unsigned &command)
{
	if (count < 2)
		throw std::runtime_error("Syntax: diamond COMMAND [OPTIONS]. To print help message: diamond help");
	const string cmd(str[1]);
	for (command = 0; command < commands_.size(); ++command)
		if (commands_[command].first == cmd || "--" + commands_[command].first == cmd)
			break;
	if (command == commands_.size())
		throw std::runtime_error("Invalid command: " + cmd + ". To print help message: diamond help");

	for (std::map<string, Option_base*>::const_iterator i = map_.begin(); i != map_.end(); ++i)
		i->second->set_default();

	vector<string> v;
	for (int i = 2; i < count; ++i) {
		if (str[i][0] == '-') {
			store_option(v);
			v.clear();
		}
		v.push_back(str[i]);
	}
	store_option(v);
}

void Command_line_parser::print_help()
{
	static const size_t col1_width = 23;
	cout << "Syntax: diamond COMMAND [OPTIONS]" << endl << endl;
	cout << "Commands:" << endl;
	for (vector<pair<string, string> >::const_iterator i = commands_.begin(); i != commands_.end(); ++i)
		if(i->second != "")
			cout << i->first << '\t' << i->second << endl;
	cout << endl;
	for (vector<const Options_group*>::const_iterator i = groups_.begin(); i != groups_.end(); ++i) {
		if ((*i)->title == "")
			continue;
		cout << (*i)->title << ":" << endl;
		for (vector<Option_base*>::const_iterator j = (*i)->options.begin(); j != (*i)->options.end(); ++j) {
			string col1 = "--" + (*j)->id;
			if ((*j)->short_id)
				col1 += string(" (-") + (*j)->short_id + ")";
			col1.append(std::max(col1_width, col1.length()) - col1.length(), ' ');
			cout << col1 << (*j)->desc << endl;
		}
		cout << endl;
	}
}
