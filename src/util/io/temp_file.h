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

#pragma once
#include "output_file.h"

struct TempFileHandler {
	void init(const char* path);
private:
	const std::string path_;
};

extern TempFileHandler temp_file_handler;

struct TempFile : public OutputFile
{

	TempFile(bool unlink = true);
	TempFile(const std::string & file_name);
	TempFile(const TempFileData& d);
	virtual void finalize() override {}
	static std::string get_temp_dir();
	static unsigned n;
	static uint64_t hash_key;
	bool unlinked;

	static TempFileData init(bool unlink);

private:

};

struct TmpFile {

	TmpFile():
		file_(nullptr)
	{
		TempFileData d = TempFile::init(true);
#ifdef _MSC_VER
		file_ = fopen(d.name.c_str(), "w+b");
#else
		file_ = fdopen(d.fd, "w+b");
#endif
		if (file_ == 0) {
			perror(("\nError opening temporary file " + d.name).c_str());
			throw std::runtime_error("Error opening temporary file " + d.name);
		}
		if (setvbuf(file_, nullptr, _IOFBF, 64 * MEGABYTES) != 0)
			throw std::runtime_error("Error calling setvbuf.");
		unlinked_ = d.unlinked;
		file_name_ = d.name;
	}

	TmpFile(const TmpFile&) = delete;
	TmpFile& operator= (const TmpFile&) = delete;
	TmpFile(TmpFile&& f) noexcept {
		file_ = f.file_;
		unlinked_ = f.unlinked_;
		file_name_ = f.file_name_;
		f.file_ = nullptr;
	}
	TmpFile& operator=(TmpFile&& f) noexcept {
		file_ = f.file_;
		unlinked_ = f.unlinked_;
		file_name_ = f.file_name_;
		f.file_ = nullptr;
		return *this;
	}

	void close() {
		if (file_) {
			fclose(file_);
			if (!unlinked_ && remove(file_name_.c_str()) != 0)
				fprintf(stderr, "Warning: Failed to delete temporary file %s\n", file_name_.c_str());
		}
		file_ = nullptr;
	}

	~TmpFile() {
		close();
	}

	void write(const void* ptr, size_t n) {
		if (fwrite(ptr, 1, n, file_) != n) {
			perror(0);
			throw std::runtime_error("Error writing to temporary file " + file_name_);
		}
	}

	void seek(int64_t p, int origin)
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

	int64_t tell()
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

	int64_t size() {
		const int64_t pos = tell();
		seek(0l, SEEK_END);
		const int64_t s = tell();
		seek(pos, SEEK_SET);
		return s;
	}

	FILE* file() {
		return file_;
	}

private:

	FILE* file_;
	bool unlinked_;
	std::string file_name_;

};