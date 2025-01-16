#include "tsv.h"
#include "../io/text_input_file.h"

using std::vector;
using std::string;
using std::numeric_limits;
using std::runtime_error;

namespace Util { namespace Tsv {

int64_t count_lines(const std::string& file_name) {
	TextInputFile f(file_name);
	int64_t n = 0;
	while (f.getline(), !f.line.empty() || !f.eof())
		++n;
	f.close();
	return n;
}

}}