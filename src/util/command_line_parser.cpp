/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
