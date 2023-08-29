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

#include <stdio.h>
#ifdef _MSC_VER
#define NOMINMAX
#include <Windows.h>
#else
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/time.h>
#endif

#include <vector>
#include "temp_file.h"
#include "../../basic/config.h"
#include "../util.h"
#include "input_file.h"

using std::vector;
using std::string;
using std::pair;

unsigned TempFile::n = 0;
uint64_t TempFile::hash_key;

#ifdef _MSC_VER
string TempFile::init(bool unlink)
#else
pair<string, int> TempFile::init(bool unlink)
#endif
{
	vector<char> buf(config.tmpdir.length() + 64);
	char *s = buf.data();
	const string prefix = config.tmpdir != "" ? config.tmpdir + dir_separator : "";

#ifdef _MSC_VER
	if (n == 0) {
		LARGE_INTEGER count;
		QueryPerformanceCounter(&count);
		hash_key = (uint64_t)(count.HighPart + count.LowPart + count.QuadPart + GetCurrentProcessId());
	}
	snprintf(s, buf.size() - 1, "%sdiamond-%llx-%u.tmp", prefix.c_str(), hash_key, n++);
#else
	snprintf(s, buf.size() - 1, "%sdiamond-tmp-XXXXXX", prefix.c_str());
#endif

#ifdef _MSC_VER
	unlinked = false;
	return string(s);
#else
	int fd = mkstemp(s);
	if (fd < 0) {
		perror(0);
		throw std::runtime_error(string("Error opening temporary file ") + string(s));
	}
	if (config.no_unlink || !unlink)
		unlinked = false;
	else
		unlinked = (::unlink(s) >= 0);
	return std::make_pair(string(s), fd);
#endif
}

TempFile::TempFile(bool unlink):
#ifdef _MSC_VER
	OutputFile(init(unlink), Compressor::NONE, "w+b")
#else
	OutputFile(init(unlink), "w+b")
#endif
{
}

TempFile::TempFile(const std::string & file_name):
	OutputFile(file_name)
{
}

string TempFile::get_temp_dir()
{
	TempFile t;
	InputFile f(t);
	f.close_and_delete();
	return extract_dir(f.file_name);
}