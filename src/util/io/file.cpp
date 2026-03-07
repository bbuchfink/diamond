/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

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