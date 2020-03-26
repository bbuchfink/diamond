// various tools related to the multiprocessing parallelization, to be moved elsewhere

#include <string>
#include <vector>
#include <sstream>
#include "multiprocessing.h"

using namespace std;

vector<string> split(const string & str, const char delim) {
	auto strstr = stringstream(str);
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