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

#include <iostream>
#include <algorithm>
#include <string.h>
#include "command_line_parser.h"

using namespace std;

OptionsGroup& CommandLineParser::add_group(const char *title, const std::vector<unsigned>& commands, bool disabled) {
	groups_.emplace_back(title, commands, disabled, this);
	return groups_.back();
}

CommandLineParser & CommandLineParser::add_command(const char * s, const char *desc, unsigned code)
{
	command_codes_[s] = code;
	commands_.emplace_back(s, desc);
	return *this;
}

void CommandLineParser::store_option(const vector<string> &v, unsigned command)
{
	if (v.size() == 0)
		return;
	OptionBase* o = nullptr;
	string id;
	vector<string> v2;

	if (v[0].length() <= 1)
		throw runtime_error("Invalid option syntax.");
	else if (v[0].substr(0, 2) == "--") {
		id = v[0].substr(2);
		map<string, OptionBase*>::const_iterator i = map_.find(id);
		if (i != map_.end())
			o = i->second;
	} else if (v[0][0] == '-') {
		id = v[0][1];
		map<char, OptionBase*>::const_iterator i = map_short_.find(v[0][1]);
		if (i != map_short_.end())
			o = i->second;
		if (v[0].length() > 2)
			v2.push_back(v[0].substr(2));
	} else
		throw runtime_error("Command line options must begin with - or --.");
	
	if (o == nullptr || o->disabled)
		throw runtime_error("Invalid option: " + id);
#ifndef EXTRA
	else if (find(o->group->commands.begin(), o->group->commands.end(), command) == o->group->commands.end())
		throw runtime_error("Option is not permitted for this workflow: " + id);
#endif
	else {
		v2.insert(v2.end(), v.begin() + 1, v.end());
		o->read(v2);
	}
}

void CommandLineParser::store(int count, const char ** str, unsigned &command)
{
	if (count < 2)
		throw runtime_error("Syntax: diamond COMMAND [OPTIONS]. To print help message: diamond help");
	string cmd(str[1]);
	if (cmd.compare(0, 2, "--") == 0)
		cmd = cmd.substr(2);
	auto it = command_codes_.find(cmd);
	if (it == command_codes_.end())
		throw runtime_error("Invalid command: " + cmd + ". To print help message: diamond help");
	command = it->second;

	for (map<string, OptionBase*>::const_iterator i = map_.begin(); i != map_.end(); ++i)
		i->second->set_default();

	vector<string> v;
	for (int i = 2; i < count; ++i) {
		if (str[i][0] == '-' && strlen(str[i]) > 1 && !isdigit(str[i][1])) {
			store_option(v, command);
			v.clear();
		}
		v.push_back(str[i]);
	}
	store_option(v, command);
}

void CommandLineParser::print_help()
{
	static const size_t col1_width = 25;
	cout << "Syntax: diamond COMMAND [OPTIONS]" << endl << endl;
	cout << "Commands:" << endl;
	for (vector<pair<string, string> >::const_iterator i = commands_.begin(); i != commands_.end(); ++i)
		if (i->second != "") {
			cout << i->first;
			for (size_t j = 0; j < col1_width - i->first.length(); ++j)
				cout << ' ';
			cout << i->second << endl;
		}
    cout << endl;
    cout << "Possible [OPTIONS] for COMMAND can be seen with syntax: diamond COMMAND" <<endl;
	cout << endl;
	cout << "Online documentation at http://www.diamondsearch.org" << endl;
}
void CommandLineParser::print_documentation(int command)
{
    static const size_t col1_width = 25;
	cout << "Options:" << endl;
    for(auto const& opt : groups_){
        for( int i = 0; i < opt.commands.size();i++) {
            if(opt.commands[i] == command){
				for (auto const& i : opt.options) {
					if (i->desc.empty())
						continue;
					string col1 = "--" + (*i).id;
					col1.append(max(col1_width, col1.length()) - col1.length(), ' ');
					cout << col1 << (*i).desc << endl;
				}
            }
        }
    }
	cout << endl;
}


void CommandLineParser::require(const char* option) {
	auto it = map_.find(option);
	if (it == map_.end())
		throw std::runtime_error("Unknown option.");
	if(!it->second->present())
		throw std::runtime_error("Missing parameter: --" + it->second->id + "/-" + it->second->short_id);
}