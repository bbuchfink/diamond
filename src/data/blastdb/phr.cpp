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

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <sstream>
#include "volume.h"
#include "util/io/mmap.h"
#include "ber.h"
#include "asn1.h"

using std::vector;
using std::string;
using std::runtime_error;
using std::out_of_range;

static std::string TagNameFromNumber(uint32_t num)
{
    static const std::unordered_map<uint32_t, string> kNames = {
        {0, "local"},        {1, "gibbsq"},      {2, "gibbmt"},    {3, "giim"},
        {4, "genbank"},      {5, "embl"},        {6, "pir"},       {7, "swissprot"},
        {8, "patent"},       {9, "other"},       {10, "general"},  {11, "gi"},
        {12, "ddbj"},        {13, "prf"},        {14, "pdb"},      {15, "tpg"},
        {16, "tpe"},        {17, "tpd"},        {18, "gpipe"},    {19, "named-annot-track"}
    };
    auto it = kNames.find(num);
    return it != kNames.end() ? it->second : ("unknown-" + std::to_string(num));
}

uint32_t Volume::read_offset(const ByteView& view, size_t base_offset, uint32_t index) const
{
    const size_t position = base_offset + static_cast<size_t>(index) * sizeof(uint32_t);
    if (position + sizeof(uint32_t) > index_.pin_length || position + sizeof(uint32_t) > view.size()) {
        throw runtime_error("Offset array exceeds PIN file size");
    }
    size_t offset = position;
    return ReadBE32(view, offset);
}

Volume::Volume(PinIndex index, MappedFile pin_mapping, MappedFile phr_mapping, MappedFile psq_mapping, int idx, OId begin, OId end)
    : pin_mapping_(std::move(pin_mapping)),
    index_(std::move(index)),
    phr_mapping_(std::move(phr_mapping)),
    psq_mapping_(std::move(psq_mapping)),
	idx(idx),
    begin(begin),
	end(end)
{
}

Volume Volume::Open(const string& path, int idx, OId begin, OId end)
{
    MappedFile pin_mapping(path + ".pin");
    PinIndex index = Volume::ParsePinFile(pin_mapping);
    MappedFile phr_mapping(path + ".phr");
    MappedFile psq_mapping(path + ".psq");
    return Volume(std::move(index), std::move(pin_mapping), std::move(phr_mapping), std::move(psq_mapping), idx, begin, end);
}

static SeqId decode_seqid(const Node& node) {
	SeqId seqid;
    for (const Node& n1 : node.children) {
        switch (n1.tag.tagNumber) {
        case 16:
            for (const Node& n2 : n1.children) {
                switch (n2.tag.tagNumber) {
                case 4:
                case 5:
                case 7:
                case 9:
                case 12:
                case 15:
                case 16:
                    for (const Node& n3 : n2.children) {
                        switch (n3.tag.tagNumber) {
                        case 16:
                            for (const Node& n4 : n3.children) {
                                switch (n4.tag.tagNumber) {
                                case 1:
                                    for (const Node& n5 : n4.children) {
                                        switch (n5.tag.tagNumber) {
                                        case 26:
                                            seqid.value.assign(n5.value.begin(), n5.value.end());
                                            break;
                                        default:
                                            ;
                                        }
                                    }
                                    break;
                                case 3:
                                    for (const Node& n5 : n4.children) {
                                        switch (n5.tag.tagNumber) {
                                        case 2:
                                            seqid.version.emplace(decode_integer(n5.value));
                                            break;
                                        default:
                                            ;
                                        }
                                    }
                                    break;
                                default:
                                    ;
                                }
                            }
                        default:
                            ;
                        }
                    }
                    break;
                case 14:
                    for (const Node& n3 : n2.children) {
                        switch (n3.tag.tagNumber) {
                        case 16:
                            for (const Node& n4 : n3.children) {
                                switch (n4.tag.tagNumber) {
                                case 0:
                                    for (const Node& n5 : n4.children) {
                                        switch (n5.tag.tagNumber) {
                                        case 26:
                                            seqid.value.assign(n5.value.begin(), n5.value.end());
                                            break;
                                        default:
                                            ;
                                        }
                                    }
                                    break;                                
                                case 3:
                                    for (const Node& n5 : n4.children) {
                                        switch (n5.tag.tagNumber) {
                                        case 26:
											seqid.chain.emplace(n5.value.begin(), n5.value.end());
                                            break;
                                        default:
                                            ;
                                        }
                                    }
                                    break;
                                default:
                                    ;
                                }
                            }
                        default:
                            ;
                        }
                    }
                default:
                    ;
                }
            }
        default:
            ;
        }
    }
    return seqid;
}

static BlastDefLine decode_defline(const Node& node, bool full_titles, bool taxids) {
    BlastDefLine defline;
    for (const Node& n1 : node.children) {
        switch (n1.tag.tagNumber) {
        case 0:
            if (full_titles) {
                for (const Node& n2 : n1.children) {
                    if (n2.tag.tagNumber == 26) {
                        defline.title.assign(n2.value.begin(), n2.value.end());
                    }
                }
            }
            break;
        case 1:
            defline.seqids.push_back(decode_seqid(n1));
            break;
        case 2:
            if (taxids) {
                for (const Node& n2 : n1.children) {
                    if (n2.tag.tagNumber == 2) {
                        defline.taxid = decode_integer(n2.value);
                    }
                }
            }
            break;
        default:
            ;
        }
    }
    return defline;
}

//#include <iostream>

static vector<BlastDefLine> decode_deflines(const char* header_data, size_t len, bool all, bool full_titles, bool taxids) {
    vector<BlastDefLine> out;
	vector<Node> nodes = decode(header_data, len);
    if (nodes.empty())
        return out;
    for (const auto& i : nodes.front().children) {
        //print_node(i, std::cout, 0);
        out.push_back(decode_defline(i, full_titles, taxids));
        if (!all && !taxids)
            break;
    }
    return out;
}

vector<BlastDefLine> Volume::deflines(std::uint32_t oid, bool all, bool full_titles, bool taxids) const
{
    if (oid >= index_.num_oids) {
        throw out_of_range("OID exceeds number of sequences in volume");
    }
    const ByteView pin_view = pin_mapping_.view();
    const size_t offset_index = static_cast<size_t>(oid);
    const ByteView phr_view = phr_mapping_.view();
    size_t header_offset = read_offset(pin_view, index_.header_offsets_offset, offset_index);
    size_t next_header_offset = read_offset(pin_view, index_.header_offsets_offset, offset_index + 1);
    if (next_header_offset < header_offset || next_header_offset > phr_view.size()) {
        throw out_of_range("Header offsets exceed PHR file size");
    }
    const size_t header_length = next_header_offset - header_offset;
	return decode_deflines(reinterpret_cast<const char*>(phr_view.data() + header_offset), header_length, all, full_titles, taxids);
 }

size_t Volume::id_len(uint32_t oid) const {
    const ByteView pin_view = pin_mapping_.view();
    const size_t offset_index = static_cast<size_t>(oid);
    return read_offset(pin_view, index_.header_offsets_offset, offset_index + 1) - read_offset(pin_view, index_.header_offsets_offset, offset_index);
}

string format_seqid(const SeqId& id) {
    if(id.value.empty()) {
        return "N/A";
	}
    std::ostringstream os;
    os << id.value;
    if (id.version) {
        os << '.' << *id.version;
    }
    if (id.chain && !id.chain->empty()) {
        os << '_' << *id.chain;
    }
    return os.str();
}

string build_title(const vector<BlastDefLine>& deflines, const char* delimiter, bool all) {
    string h;
    for (auto i = deflines.cbegin(); i != deflines.cend(); ++i) {
        if (i != deflines.cbegin()) {
            if (!all)
                break;
            h += delimiter;
        }
        h += format_seqid(i->seqids.front());
        h += ' ';
        h += i->title;
        
    }
    if (h.empty())
        h = "N/A";
    return h;
}