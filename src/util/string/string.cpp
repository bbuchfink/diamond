#include <sstream>
#include <iomanip>
#include "string.h"
#include "fixed_string.h"

using std::string;
using std::to_string;
using std::stringstream;
using std::setprecision;
using std::fixed;

std::array<char, 16> fixed_string_seed;

string convert_size(size_t size) {
	static const char *SIZES[] = { "B", "KB", "MB", "GB" };
	size_t div = 0;
	size_t rem = 0;

	while (size >= 1024 && div < (sizeof SIZES / sizeof *SIZES)) {
		rem = (size % 1024);
		div++;
		size /= 1024;
	}

	std::stringstream ss;
	ss << std::fixed << std::setprecision(1) << (double)size + (double)rem / 1024.0 << ' ' << SIZES[div];
	return ss.str();
}

namespace Util { namespace String {

string replace(const std::string& s, char a, char b) {
	string r = s;
	size_t i = r.find_first_of(a);
	while (i != string::npos) {
		r[i++] = b;
		i = r.find_first_of(a, i);
	}
	return r;
}

std::string ratio_percentage(const double x, const double y) {
	stringstream ss;
	ss << fixed << setprecision(0) << x << '/' << y << " (" << setprecision(2) << x / y * 100.0 << "%)";
	return ss.str();
}

std::string ratio_percentage(const size_t x, const size_t y) {
	return ratio_percentage((double)x, (double)y);
}

int64_t interpret_number(const std::string& s) {
	stringstream ss(s);
	double n;
	ss >> n;
	char c = 0;
	ss >> c;
	if (ss.eof())
		throw std::runtime_error("Missing size specifier in number: " + s + ". Permitted values: 'G'");
	if (c != 'G')
		throw std::runtime_error(string("Invalid size specifier (") + c + ") in number: " + s + ". Permitted values: 'G'");
	ss >> c;
	if (!ss.eof())
		throw std::runtime_error("Invalid number format: " + s);
	return int64_t(n * 1e9);
}

}}