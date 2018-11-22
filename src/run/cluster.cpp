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

#include <string>
#include "../util/system/system.h"
#include "../util/util.h"

using namespace std;

void build_homology_graph(const string &input_file) {
	FILE * stream;
	stream = popen(join(' ', { executable_path(), "blastp",
		"-q", input_file,
		-d, "r");

	if (stream) {
		while (!feof(stream))
			if (fgets(buffer, max_buffer, stream) != nullptr)
				out.append(buffer);
		if (pclose(stream))
			throw runtime_error("Error running command " + cmd);
	}
	else
		throw runtime_error("Error running command " + cmd);
}

void cluster() {

}