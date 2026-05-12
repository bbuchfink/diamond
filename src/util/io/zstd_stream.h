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
#include <zstd.h>
#include "stream_entity.h"
#include "input_stream_buffer.h"

struct ZstdSink : public StreamEntity
{
	ZstdSink(StreamEntity* prev);
	virtual void close() override;
	virtual void write(const char* ptr, size_t count) override;
private:
	ZSTD_CStream* stream;
};

struct ZstdSource : public StreamEntity
{
	ZstdSource(InputStreamBuffer* prev);
	virtual size_t read(char* ptr, size_t count) override;
	virtual void close() override;
	virtual void rewind() override;
	virtual bool eof() override;
private:
	void init();
	ZSTD_DStream* stream;
	ZSTD_inBuffer in_buf;
	bool eos_;
};

size_t zstd_decompress(FILE* src, void* dst, size_t dstCapacity);