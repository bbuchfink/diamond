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
#include <optional>
#include <string>
#include <vector>
#include "util/io/mmap.h"
#include "basic/value.h"

struct SeqId {
    std::string type;
    std::string value;
    std::optional<std::int64_t> version;
    std::optional<std::string> chain;
};

struct BlastDefLine {
    std::string title;
    std::vector<SeqId> seqids;
    std::optional<std::int64_t> taxid;
};

struct PinIndex {
    std::uint32_t version = 0;
    bool is_protein = false;
    std::uint32_t volume_number = 0; // only meaningful for version 5
    std::string title;
    std::string lmdb_file; // version 5 only
    std::string date;
    std::uint32_t num_oids = 0;
    std::uint64_t total_length = 0;
    std::uint32_t max_length = 0;
    std::size_t header_offsets_offset = 0;
    std::size_t sequence_offsets_offset = 0;
    std::size_t ambiguity_offsets_offset = 0; // nucleotide only
    std::size_t pin_length = 0;
};

std::string build_title(const std::vector<BlastDefLine>& deflines, const char* delimiter, bool all);
std::string format_seqid(const SeqId& id);

class Volume {
public:
    Volume() = default;
    static Volume Open(const std::string& path, int idx, OId begin, OId end);

    const PinIndex& index() const { return index_; }

    std::vector<BlastDefLine> deflines(uint32_t oid, bool all, bool full_titles, bool taxids) const;
    std::vector<Letter> sequence(uint32_t oid) const;
	Loc length(uint32_t oid) const;
    size_t id_len(uint32_t oid) const;
    int idx;
    OId begin = 0, end = 0;

private:
    Volume(PinIndex index, MappedFile pin_mapping, MappedFile phr_mapping, MappedFile psq_mapping, int idx, OId begin, OId end);
    static PinIndex ParsePinFile(const MappedFile& mapping);
    uint32_t read_offset(const MappedFile::View& view, size_t base_offset, uint32_t index) const;

    MappedFile pin_mapping_;
    PinIndex index_;
    MappedFile phr_mapping_;
    MappedFile psq_mapping_;
    
};