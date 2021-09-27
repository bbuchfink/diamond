#pragma once

#include <string>
#include "../io/text_input_file.h"

namespace Util { namespace Tsv {

std::string fetch_block(TextInputFile& f, std::string& buf);
std::string column(const std::string& line, const size_t i);
std::string columns(const std::string& line, const size_t begin, const size_t end);
size_t column_count(const std::string& line);
std::vector<std::string> extract_column(const std::string& buf, const size_t i);

}}