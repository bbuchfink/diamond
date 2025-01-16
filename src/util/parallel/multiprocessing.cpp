/****
DIAMOND protein aligner
Copyright (C) 2019-2024 Max Planck Society for the Advancement of Science e.V.

Code developed by Klaus Reuter

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

// various tools related to the multiprocessing parallelization, to be moved elsewhere

#include <string>
#include <vector>
#include <fstream>
#include "multiprocessing.h"

using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::string;
using std::vector;

vector<string> split(const string & str, const char delim) {
	stringstream strstr(str);
	string segment;
	vector<string> segments;
	while(getline(strstr, segment, delim))
	{
		segments.push_back(segment);
	}
	return segments;
}

string join(const vector<string> & tokens, const char delim) {
	string buf;
	for (size_t i=0; i<tokens.size(); i++) {
		buf.append(tokens[i]);
		if (i+1 < tokens.size()) {
			buf += delim;
		}
	}
	return buf;
}

string quote(const string & str) {
	return '"' + str + '"';
}

string unquote(const string & str) {
	if ((str.size() >= 2) && (str.front() == '"') && (str.back() == '"')) {
		return string(str.begin()+1, str.end()-1);
	} else {
		return str;
	}
}

void copy(const string & src_file_name, const string & dst_file_name) {
	ifstream src(src_file_name, ios::binary);
	ofstream dst(dst_file_name, ios::binary);
	dst << src.rdbuf();
}

string join_path(const string & path_1, const string & path_2) {
#ifdef __unix
    const string sep = "/";
#else
    const string sep = "\\";
#endif
    return path_1 + sep + path_2;
}

bool file_exists(const string & file_name) {
    ifstream fp(file_name.c_str());
    return fp.good();
}
