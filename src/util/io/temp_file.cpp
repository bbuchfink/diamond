/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#include <stdio.h>
#ifdef _MSC_VER
#define NOMINMAX
#define NT_SUCCESS(Status) (((NTSTATUS)(Status)) >= 0)
#define STATUS_UNSUCCESSFUL ((NTSTATUS)0xC0000001L)
#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <bcrypt.h>
#else
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/time.h>
#endif

#include <iostream>
#include <array>
#include <vector>
#include "temp_file.h"
#include "basic/config.h"
#include "../util.h"
#include "input_file.h"

using std::vector;
using std::string;
using std::pair;
using std::array;
using std::runtime_error;
using std::to_string;

TempFileHandler temp_file_handler;

void TempFileHandler::init(const char* path) {
	if (!path_.empty())
		throw runtime_error("Double init of TempFileHandler");
#ifdef _MSC_VER
	NTSTATUS status = STATUS_UNSUCCESSFUL;
	array<UCHAR, 20> buf;
	
	if(!NT_SUCCESS(status = BCryptGenRandom(nullptr, buf.data(), 20, BCRYPT_USE_SYSTEM_PREFERRED_RNG)))
		throw runtime_error("Error " + to_string(status) + " returned by BCryptGenRandom");
#endif
}

unsigned TempFile::n = 0;
uint64_t TempFile::hash_key;

TempFileData TempFile::init(bool unlink)
{
	TempFileData r;
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
	r.unlinked = false;
	r.name = s;
	return r;
#else
	int fd = mkstemp(s);
	if (fd < 0) {
		perror(0);
		throw std::runtime_error(string("Error opening temporary file ") + string(s));
	}
	if (config.no_unlink || !unlink)
		r.unlinked = false;
	else
		r.unlinked = (::unlink(s) >= 0);
	r.name = s;
	r.fd = fd;
	return r;
#endif
}

TempFile::TempFile(const TempFileData& d):
	OutputFile(d, Compressor::NONE, "w+b"),
	unlinked(d.unlinked)
{
}

TempFile::TempFile(bool unlink):
	TempFile(init(unlink))
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