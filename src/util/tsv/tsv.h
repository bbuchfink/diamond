#pragma once

#include <string>
#include "file.h"
#include "file_read.h"

namespace Util { namespace Tsv {

int64_t count_lines(const std::string& file_name);

}}