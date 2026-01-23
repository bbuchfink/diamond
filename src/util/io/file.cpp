/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#include <errno.h>
#include <string.h> // strerror
#include "file.h"
#include "temp_file.h"
#include "util/data_structures/mem_buffer.h"

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

File::File(const string& name) {
	file_ = fopen(name.c_str(), "r+b");
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
		if (auto_delete_ && !unlinked_ && remove(file_name_.c_str()) != 0)
			fprintf(stderr, "Warning: Failed to delete temporary file %s\n", file_name_.c_str());
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