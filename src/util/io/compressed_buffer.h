/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

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