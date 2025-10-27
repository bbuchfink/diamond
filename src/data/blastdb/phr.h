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
#include <string>
#include <vector>
#include <stdexcept>
#include "util/io/mmap.h"
#include "basic/value.h"

// ---- replace your PinIndex with this ----
struct PinIndex {
    uint32_t version = 0;            // 4 or 5
    uint32_t dbtype = 0;            // 0: DNA, 1: Protein
    std::string_view title;
    std::string_view timestamp;
    uint32_t nseq = 0;

    uint64_t residue_count = 0;
    uint32_t max_seq_len = 0;

    std::vector<uint64_t> hdr_offsets;
};

struct Phr {

    Phr(const std::string& phr, const std::string& pin);

    bool parse_record(OId i,
        std::vector<std::pair<const uint8_t*, size_t>>& titles,
        std::vector<std::string>& ids,
        std::vector<uint64_t>& taxids);

    size_t size() const {
        return idx_.nseq;
    }

    size_t len(OId i) const {
        if (i + 1 >= (OId)idx_.hdr_offsets.size())
            throw std::out_of_range("Phr::len");
        return (size_t)(idx_.hdr_offsets[i + 1] - idx_.hdr_offsets[i]);
	}

    ~Phr();

private:

    PinIndex idx_;
    MappedFile phr_;

};