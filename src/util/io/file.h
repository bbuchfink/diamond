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

#pragma once
#include <cstdint>
#include <memory>
#include <string>
#include <cstdio>
#include <functional>
#include <string.h>
#include "util/enum.h"
#include "decompressor.h"
#include "compressor.h"
#include "util/memory/memory_resource.h"

const size_t MEGABYTES = 1 << 20;
const size_t GIGABYTES = 1 << 30;
const size_t KILOBYTES = 1 << 10;

struct Temporary {};

struct File {

	enum struct Flags : int { DETECT_COMPRESSION = 1, TREAT_BLANK_AS_STDIN = 2, TREAT_BLANK_AS_STDOUT = 4, NONE = 0 };

	File(Temporary);
	File(const std::string& name, const char* mode, Flags flags = Flags::NONE, CompressionLib compression = CompressionLib::NONE);
	File(const File&) = delete;
	File& operator= (const File&) = delete;
	File(File&& f) noexcept;
	File& operator=(File&& f) noexcept;
	void close();
	~File();
	void seek(int64_t p, int origin);
	void seek(int64_t p);
	void rewind();
	int64_t tell();
	int64_t size();
	FILE* file();
	const std::string& name() const {
		return file_name_;
	}
	void read(void* ptr, size_t n);
	size_t read_max(void* ptr, size_t n);
	const char* read_bytes(size_t n);
	template<typename T>
	void read(T& x) {
		read(&x, sizeof(T));
	}
	void read_c_str(std::string& s);
	void read_c_str(std::pmr::string& s);

	void write(const void* ptr, size_t n);
	template<typename T>
	void write(const T& x) {
		write(&x, sizeof(T));
	}
	void write_c_str(const char* s);
	
	bool eof() const;
	std::string peek(int64_t n);
	const char* getline();
	const char* getdelim(char delimiter);
	size_t read_raw(std::string& dst, size_t count);

	template<typename It>
	void read_to(It dst, char delimiter) {
		const ssize_t n = decompressor_->getdelim(&line_buf_, &line_buf_size_, delimiter, file_);
		if (ferror(file_))
			throw std::runtime_error("Error reading file " + file_name_ + ". " + strerror(errno));
		std::copy(line_buf_, line_buf_ + n, dst);
	}

	uint64_t read_to_fasta_record_end(std::string& dst);
	static int64_t count_lines(const std::string& file_name);
	void read_text_mt(int64_t max_size, int threads, std::function<void(int64_t chunk, const char*, const char*)>& callback);
	size_t line_count = 0;

private:

	FILE* file_ = nullptr;
	bool auto_delete_, unlinked_, seekable_;
	std::string file_name_;
	char* line_buf_;
	size_t line_buf_size_;
	std::unique_ptr<Decompressor> decompressor_;
	std::unique_ptr<CompressorX> compressor_;

};

DEFINE_ENUM_OPERATORS(File::Flags)