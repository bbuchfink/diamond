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
#include "util/io/file.h"

inline uint32_t ReadBE32(File& file)
{
    uint8_t buf[4];
    file.read(buf, 4);
    uint32_t value = 0;
    for (int i = 0; i < 4; ++i) {
        value = (value << 8) | buf[i];
    }
    return value;
}

inline uint64_t ReadLE64(File& file)
{
    uint8_t buf[8];
    file.read(buf, 8);
    uint64_t value = 0;
    for (int i = 7; i >= 0; --i) {
        value = (value << 8) | static_cast<uint64_t>(buf[i]);
    }
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

inline std::string ReadPascalString(File& file) {
	const uint32_t l = ReadBE32(file);
    std::string s(l, '\0');
    file.read(&s.front(), l);
    return s;
}