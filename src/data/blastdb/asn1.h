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
#include <stdexcept>
#include <vector>

enum class Class : uint8_t {
    Universal = 0,
    Application = 1,
    ContextSpecific = 2,
    Private = 3,
};

enum class UniversalTag : uint8_t {
    Eoc = 0,
    Boolean = 1,
    Integer = 2,
    BitString = 3,
    OctetString = 4,
    Null = 5,
    ObjectIdentifier = 6,
    Utf8String = 12,
    Sequence = 16,
    Set = 17,
    PrintableString = 19,
    T61String = 20,
    Ia5String = 22,
    UtcTime = 23,
    GeneralizedTime = 24,
    BmpString = 30,
};

struct TagInfo {
    Class tagClass{};
    bool constructed{};
    uint32_t tagNumber{};
};

struct Node {
    TagInfo tag;
    std::vector<uint8_t> value;
    std::vector<Node> children;
};

class DecodeError : public std::runtime_error {
  public:
    using std::runtime_error::runtime_error;
};

std::vector<Node> decode(const char *data, std::size_t length);
void print_node(const Node& node, std::ostream& os, int depth = 0);