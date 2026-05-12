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
#include <vector>

struct z_stream_s;
struct ZSTD_CCtx_s;

struct CompressedBuffer {

	static constexpr size_t BUF_SIZE = 32768;

	CompressedBuffer();
	~CompressedBuffer();
	void write(const char* ptr, size_t n);

	template<typename T>
	void write(const T& x) {
		write(reinterpret_cast<const char*>(&x), sizeof(T));
	}

	void finish();
	void clear();

	const char* data() const {
		return buf_.data();
	}

	size_t size() const {
		return size_;
	}

private:
	std::vector<char> buf_;
#ifdef WITH_ZSTD
	ZSTD_CCtx_s* stream_;
#else
	z_stream_s* stream_;
#endif
	size_t size_;

};