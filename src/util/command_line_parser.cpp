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

#include <algorithm>
#include <string.h>
#include <iostream>
#include "command_line_parser.h"

using std::string;
using std::vector;
using std::map;
using std::runtime_error;
using std::cout;
using std::max;
using std::endl;
using std::pair;

OptionsGroup& CommandLineParser::add_group(const char *title, const std::vector<int>& commands, bool disabled) {
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
		for (size_t i = 0; i < opt.commands.size(); i++) {
			if (opt.commands[i] == command) {
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