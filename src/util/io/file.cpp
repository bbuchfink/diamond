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

#include <errno.h>
#include <string.h> // strerror
#include "file.h"
#include "temp_file.h"
#include "util/data_structures/mem_buffer.h"
#include "util/system/system.h"

using std::string;
using std::runtime_error;

File::File(Temporary) {	
	TempFileData d = TempFile::init(true);
#ifdef _MSC_VER
	file_ = fopen(d.name.c_str(), "w+b");
#else
	file_ = fdopen(d.fd, "w+b");
#endif
	if (file_ == 0)
		throw runtime_error("Error opening temporary file " + d.name + ". " + strerror(errno));
	if (setvbuf(file_, nullptr, _IOFBF, 64 * MEGABYTES) != 0)
		throw runtime_error("Error setting buffer size for temporary file " + d.name + ". " + strerror(errno));
	unlinked_ = d.unlinked;
	file_name_ = d.name;
	auto_delete_ = true;
}

File::File(const string& name, const char* mode) {
	file_ = fopen(name.c_str(), mode);
	if (file_ == 0)
		throw runtime_error("Error opening file " + name + ". " + strerror(errno));
	unlinked_ = false;
	file_name_ = name;
	auto_delete_ = false;
}

File::File(File&& f) noexcept {
	close();
	file_ = f.file_;
	unlinked_ = f.unlinked_;
	file_name_ = f.file_name_;
	auto_delete_ = f.auto_delete_;
	f.file_ = nullptr;
}

File& File::operator=(File&& f) noexcept {
	close();
	file_ = f.file_;
	unlinked_ = f.unlinked_;
	file_name_ = f.file_name_;
	auto_delete_ = f.auto_delete_;
	f.file_ = nullptr;
	return *this;
}

void File::close() {
	if (file_) {
		fclose(file_);
		if (auto_delete_ && !unlinked_)
			remove_tmp_file(file_name_);
	}
	file_ = nullptr;
}

File::~File() {
	close();
}

void File::write(const void* ptr, size_t n) {
	if (fwrite(ptr, 1, n, file_) != n) {
		perror(0);
		throw std::runtime_error("Error writing to temporary file " + file_name_);
	}
}

void File::seek(int64_t p, int origin)
{
#ifdef WIN32
	if (_fseeki64(file_, p, origin) != 0) {
		perror(0);
		throw std::runtime_error("Error calling fseek.");
	}
#else
	if (fseek(file_, p, origin) != 0) {
		perror(0);
		throw std::runtime_error("Error calling fseek.");
	}
#endif
}

int64_t File::tell()
{
#ifdef WIN32
	int64_t x;
	if ((x = _ftelli64(file_)) == (int64_t)-1)
		throw std::runtime_error("Error executing ftell on stream " + file_name_);
	return x;
#else
	const long n = ftell(file_);
	if (n < 0) {
		perror(0);
		throw std::runtime_error("Error calling ftell.");
	}
	return n;
#endif
}

int64_t File::size() {
	const int64_t pos = tell();
	seek(0l, SEEK_END);
	const int64_t s = tell();
	seek(pos, SEEK_SET);
	return s;
}

FILE* File::file() {
	return file_;
}

void File::read(void* ptr, size_t n) {
	const size_t r = fread(ptr, 1, n, file_);
	if(r != n)
		throw runtime_error("Error reading file " + file_name_ + ". " + strerror(errno));
}

size_t File::read_max(void* ptr, size_t n) {
	return fread(ptr, 1, n, file_);
}

const char* File::read(size_t n) {
	static MemBuffer<char> buf;
	buf.resize(n);
	read(buf.begin(), n);	
	return buf.begin();
}