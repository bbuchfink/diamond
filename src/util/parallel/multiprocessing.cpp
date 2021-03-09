// various tools related to the multiprocessing parallelization, to be moved elsewhere

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "multiprocessing.h"

using namespace std;

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
