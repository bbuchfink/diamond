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
#include <string>

struct Temporary {};

struct File {

	File(Temporary);
	File(const std::string& name, const char* mode);
	File(const File&) = delete;
	File& operator= (const File&) = delete;
	File(File&& f) noexcept;
	File& operator=(File&& f) noexcept;
	void close();
	~File();
	void write(const void* ptr, size_t n);
	void seek(int64_t p, int origin);
	int64_t tell();
	int64_t size();
	FILE* file();
	void read(void* ptr, size_t n);
	size_t read_max(void* ptr, size_t n);
	const char* read(size_t n);
	template<typename T>
	void read(T& x) {
		read(&x, sizeof(T));
	}
	template<typename T>
	void write(const T& x) {
		write(&x, sizeof(T));
	}

private:

	FILE* file_ = nullptr;
	bool auto_delete_, unlinked_;
	std::string file_name_;

};