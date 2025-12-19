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
#include "util/io/mmap.h"

using Byte = std::uint8_t;
using ByteView = MappedFile::View;

inline std::uint32_t ReadBE32(const ByteView& buffer, std::size_t& offset)
{
    if (offset + 4 > buffer.size()) {
        throw std::runtime_error("Unexpected end of file while reading 32-bit value");
    }

    std::uint32_t value = 0;
    for (int i = 0; i < 4; ++i) {
        value = (value << 8) | buffer[offset + i];
    }
    offset += 4;
    return value;
}

inline static std::uint64_t ReadLE64(const ByteView& buffer, std::size_t& offset)
{
    if (offset + 8 > buffer.size()) {
        throw std::runtime_error("Unexpected end of file while reading 64-bit value");
    }

    std::uint64_t value = 0;
    for (int i = 7; i >= 0; --i) {
        value = (value << 8) | static_cast<std::uint64_t>(buffer[offset + i]);
    }
    offset += 8;
    return value;
}

inline int64_t decode_integer(const std::vector<std::uint8_t>& data) {
    if (data.size() > sizeof(int64_t)) {
        return {};
    }
    int64_t value = (data[0] & 0x80) ? -1 : 0;
    for (uint8_t byte : data) {
        value = (value << 8) | byte;
    }
    return value;
}

inline std::string ReadPascalString(const ByteView& buffer, std::size_t& offset)
{
    const std::uint32_t length = ReadBE32(buffer, offset);
    if (offset + length > buffer.size()) {
        throw std::runtime_error("String length exceeds file size");
    }

    std::string result(reinterpret_cast<const char*>(buffer.data() + offset), length);
    offset += length;
    return result;
}