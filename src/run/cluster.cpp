/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#include <stdexcept>
#include <string>
#include <stdio.h>
#include "../util/system/system.h"
#include "../util/util.h"
#include "../basic/config.h"

using namespace std;

void build_homology_graph(const string &input_file) {
	FILE * stream;
	stream = POPEN(join(' ', { executable_path(), "blastp",
		"-q", input_file,
		"-d", input_file,
		"-c1" }).c_str(), "r");
	if (stream == NULL)
		throw runtime_error("Error executing popen.");
	char buffer[512];
	while (!feof(stream) && !ferror(stream)) {
		if (fgets(buffer, 512, stream) != nullptr)
			cout << buffer;
	}
	if (PCLOSE(stream))
		throw runtime_error("Error executing pclose.");
}

void cluster() {
	build_homology_graph(config.query_file);
}