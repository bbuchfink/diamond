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

#include <algorithm>
#include <sstream>
#include <iomanip>
#include <ios>
#include "asn1.h"
#include "ber.h"

using std::vector;
using std::string;

constexpr uint8_t kClassMask = 0xc0;
constexpr uint8_t kConstructedMask = 0x20;
constexpr uint8_t kShortTagMask = 0x1f;

static TagInfo parse_tag(const char *data, size_t length, size_t &offset) {
    if (offset >= length) {
        throw DecodeError{"unexpected end of buffer while reading tag"};
    }

    uint8_t first = static_cast<uint8_t>(data[offset++]);
    TagInfo info{};
    info.tagClass = static_cast<Class>((first & kClassMask) >> 6);
    info.constructed = (first & kConstructedMask) != 0;

    uint8_t tag = first & kShortTagMask;
    if (tag != kShortTagMask) {
        info.tagNumber = tag;
        return info;
    }

    info.tagNumber = 0;
    bool continueReading = true;
    int shiftCount = 0;
    while (continueReading) {
        if (offset >= length) {
            throw DecodeError{"unexpected end of buffer while reading long tag"};
        }
        uint8_t byte = static_cast<uint8_t>(data[offset++]);
        continueReading = (byte & 0x80) != 0;
        info.tagNumber = (info.tagNumber << 7) | (byte & 0x7F);
        shiftCount += 7;

        if (shiftCount > 28) {
            throw DecodeError{"tag number is excessively large"};
        }
    }
    return info;
}

static size_t parse_length(const char *data, size_t length, size_t &offset, bool &indefinite) {
    if (offset >= length) {
        throw DecodeError{"unexpected end of buffer while reading length"};
    }
    uint8_t first = static_cast<uint8_t>(data[offset++]);
    if ((first & 0x80) == 0) {
        indefinite = false;
        return first;
    }
    uint8_t count = first & 0x7F;
    if (count == 0) {
        indefinite = true;
        return 0;
    }
    if (count > sizeof(size_t)) {
        throw DecodeError{"length uses more bytes than supported"};
    }
    if (offset + count > length) {
        throw DecodeError{"unexpected end of buffer while reading long length"};
    }
    size_t value = 0;
    for (uint8_t i = 0; i < count; ++i) {
        value = (value << 8) | static_cast<uint8_t>(data[offset++]);
    }
    indefinite = false;
    return value;
}

static bool is_eoc(const char *data, size_t length, size_t offset) {
    return offset + 1 < length && data[offset] == 0 && data[offset + 1] == 0;
}

static vector<Node> decode_impl(const char *data, size_t length, size_t &offset, bool stop_at_eoc) {
    vector<Node> nodes;

    while (offset < length) {
        if (stop_at_eoc && is_eoc(data, length, offset)) {
            offset += 2;
            break;
        }

        TagInfo tag = parse_tag(data, length, offset);

        bool indefinite = false;
        size_t contentLength = parse_length(data, length, offset, indefinite);

        if (!indefinite && offset + contentLength > length) {
            throw DecodeError{"content length exceeds available data"};
        }

        Node node{};
        node.tag = tag;

        if (tag.constructed) {
            size_t end = indefinite ? length : offset + contentLength;
            node.children = decode_impl(data, end, offset, indefinite);

            if (!indefinite && offset != end) {
                throw DecodeError{"constructed element did not consume its content"};
            }
        } else {
            size_t end = indefinite ? length : offset + contentLength;
            if (indefinite) {
                throw DecodeError{"indefinite length used for primitive value"};
            }
            node.value.insert(node.value.end(), data + offset, data + end);
            offset = end;
        }

        nodes.push_back(std::move(node));
    }

    return nodes;
}

vector<Node> decode(const char *data, size_t length) {
    if (data == nullptr && length > 0) {
        throw DecodeError{"null data pointer with non-zero length"};
    }
    std::size_t offset = 0;
    return decode_impl(data, length, offset, false);
}

static std::string class_label(Class cls) {
    switch (cls) {
    case Class::Universal:
        return "Universal";
    case Class::Application:
        return "Application";
    case Class::ContextSpecific:
        return "Context-specific";
    case Class::Private:
        return "Private";
    }
    return "Unknown";
}

static std::string universal_tag_name(uint32_t tagNumber) {
    switch (static_cast<UniversalTag>(tagNumber)) {
    case UniversalTag::Eoc:
        return "EOC";
    case UniversalTag::Boolean:
        return "BOOLEAN";
    case UniversalTag::Integer:
        return "INTEGER";
    case UniversalTag::BitString:
        return "BIT STRING";
    case UniversalTag::OctetString:
        return "OCTET STRING";
    case UniversalTag::Null:
        return "NULL";
    case UniversalTag::ObjectIdentifier:
        return "OBJECT IDENTIFIER";
    case UniversalTag::Utf8String:
        return "UTF8String";
    case UniversalTag::Sequence:
        return "SEQUENCE";
    case UniversalTag::Set:
        return "SET";
    case UniversalTag::PrintableString:
        return "PrintableString";
    case UniversalTag::T61String:
        return "T61String";
    case UniversalTag::Ia5String:
        return "IA5String";
    case UniversalTag::UtcTime:
        return "UTCTime";
    case UniversalTag::GeneralizedTime:
        return "GeneralizedTime";
    case UniversalTag::BmpString:
        return "BMPString";
    }
    return "";
}

static std::string hex_dump(const std::vector<std::uint8_t>& data) {
    std::ostringstream oss;
    for (std::uint8_t b : data) {
        oss << std::hex << std::setw(2) << std::setfill('0')
            << static_cast<int>(b) << ' ';
    }
    std::string out = oss.str();
    if (!out.empty()) {
        out.pop_back();
    }
    return out;
}

static bool is_printable_ascii(const std::vector<std::uint8_t>& data) {
    return std::all_of(data.begin(), data.end(), [](std::uint8_t b) {
        return b == 0x0A || (b >= 0x20 && b <= 0x7E);
        });
}

static std::string decode_oid(const std::vector<std::uint8_t>& data) {
    if (data.empty()) {
        return {};
    }

    std::ostringstream oss;
    uint8_t first = data.front();
    oss << static_cast<int>(first / 40) << '.' << static_cast<int>(first % 40);

    std::uint32_t value = 0;
    bool seen = false;
    for (std::size_t i = 1; i < data.size(); ++i) {
        std::uint8_t byte = data[i];
        value = (value << 7) | (byte & 0x7F);
        seen = true;
        if ((byte & 0x80) == 0) {
            oss << '.' << value;
            value = 0;
            seen = false;
        }
    }

    if (seen) {
        return {};
    }

    return oss.str();
}

static std::string decode_string(const std::vector<std::uint8_t>& data) {
    return std::string(data.begin(), data.end());
}

static std::string describe_value(const Node& node) {
    if (node.tag.constructed) {
        return {};
    }

    if (node.tag.tagClass != Class::Universal) {
        return {};
    }

    switch (static_cast<UniversalTag>(node.tag.tagNumber)) {
    case UniversalTag::Boolean:
        if (node.value.size() == 1) {
            return node.value[0] ? "TRUE" : "FALSE";
        }
        return {};
    case UniversalTag::Integer: {
        std::string decoded = std::to_string(decode_integer(node.value));
        return decoded.empty() ? std::string{} : decoded;
    }
    case UniversalTag::OctetString:
        if (is_printable_ascii(node.value)) {
            return '"' + decode_string(node.value) + '"';
        }
        return {};
    case UniversalTag::Null:
        return "NULL";
    case UniversalTag::ObjectIdentifier: {
        std::string oid = decode_oid(node.value);
        return oid.empty() ? std::string{} : oid;
    }
    case UniversalTag::Utf8String:
    case UniversalTag::PrintableString:
    case UniversalTag::T61String:
    case UniversalTag::Ia5String:
    case UniversalTag::BmpString:
        return decode_string(node.value);
    case UniversalTag::UtcTime:
    case UniversalTag::GeneralizedTime:
        return decode_string(node.value);
    case UniversalTag::BitString:
    case UniversalTag::Sequence:
    case UniversalTag::Set:
    case UniversalTag::Eoc:
        return {};
    }
    return {};
}

void print_node(const Node& node, std::ostream& os, int depth) {
    string indent(static_cast<size_t>(depth), '-');

    os << indent << "Class: " << class_label(node.tag.tagClass)
        << ", Constructed: " << std::boolalpha << node.tag.constructed
        << ", Tag: " << node.tag.tagNumber;

    if (node.tag.tagClass == Class::Universal) {
        string name = universal_tag_name(node.tag.tagNumber);
        if (!name.empty()) {
            os << " (" << name << ')';
        }
    }
    os << '\n';

    if (!node.value.empty()) {
        string decoded = describe_value(node);
        if (!decoded.empty()) {
            os << indent << "  Decoded: " << decoded << '\n';
        }
        else
            os << indent << "  String: " << string(node.value.begin(), node.value.end()) << '\n';
        os << indent << "  Raw (" << node.value.size() << " bytes): "
            << hex_dump(node.value) << '\n';
    }

    for (const auto& child : node.children) {
        print_node(child, os, depth + 1);
    }
}