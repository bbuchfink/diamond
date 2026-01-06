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