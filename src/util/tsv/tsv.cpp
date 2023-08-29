#include "tsv.h"
#include "../string/tokenizer.h"
#include "../io/text_input_file.h"

using std::vector;
using std::string;
using std::numeric_limits;
using std::runtime_error;

namespace Util { namespace Tsv {

string fetch_block(TextInputFile& f, string& buf) {
	string key, k;
	f.getline();
	if (f.line.empty())
		return string();
	String::Tokenizer tok(f.line, "\t");
	tok >> key;
	if (key.empty())
		throw std::runtime_error("Empty key in TSV file.");
	buf = f.line;
	while (f.getline(), !f.eof() || !f.line.empty()) {
		tok.set(f.line.data());
		tok >> k;
		if (k != key) {
			f.putback_line();
			return key;
		}
		buf.append("\n");
		buf.append(f.line);
	}
	return key;
}

std::string column(const std::string& line, const size_t i)
{
	String::Tokenizer tok(line, "\t");
	string s;
	for (size_t j = 0; j < i; ++j)
		tok >> String::Skip();
	tok >> s;
	return s;
}

std::string columns(const std::string& line, const size_t begin, const size_t end)
{
	String::Tokenizer tok(line, "\t");
	for (size_t j = 0; j < begin; ++j)
		tok >> String::Skip();
	for (size_t j = begin; j < end && tok.good(); ++j);
	return {};
}

size_t column_count(const std::string& line)
{
	if (line.empty())
		return 0;
	size_t n = 1, i = 0;
	while ((i = line.find_first_of('\t', i)) != string::npos) {
		++n;
		++i;
	}
	return n;
}

std::vector<std::string> extract_column(const std::string& buf, const size_t i)
{
	Util::String::Tokenizer tok(buf, "\n");
	string l;
	vector<string> out;
	while (tok.good() && (tok >> l, !l.empty())) {
		out.push_back(column(l, i));
	}
	return out;
}

int64_t count_lines(const std::string& file_name) {
	TextInputFile f(file_name);
	int64_t n = 0;
	while (f.getline(), !f.line.empty() || !f.eof())
		++n;
	f.close();
	return n;
}

template<>
int64_t convert_string<int64_t>(const char* s) {
	char* end;
	long long i = strtoll(s, &end, 10);
	if ((i == 0 && strcmp(s, "0") != 0) || i == LLONG_MAX || i == LLONG_MIN || *end != '\0' || i < numeric_limits<int64_t>::min() || i > numeric_limits<int64_t>::max())
		throw runtime_error(string("Error converting integer value: ") + s);
	return i;
}

template<>
int32_t convert_string<int32_t>(const char* s) {
	const int64_t i = convert_string<int64_t>(s);
	if(i < (int64_t)INT32_MIN || i > (int64_t)INT32_MAX)
		throw runtime_error(string("Error converting integer value: ") + s);
	return (int32_t)i;
}

}}