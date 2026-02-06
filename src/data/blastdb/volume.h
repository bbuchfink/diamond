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
#include "basic/value.h"
#include "util/io/file.h"
#include "util/optional.h"
#include "../sequence_file.h"

struct SeqId {
    std::string type;
    std::string value;
    optional<int64_t> version;
    optional<std::string> chain;
};

struct BlastDefLine {
    std::string title;
    std::vector<SeqId> seqids;
    optional<TaxId> taxid;
};

struct PinIndex {
    uint32_t version = 0;
    bool is_protein = false;
    uint32_t volume_number = 0; // only meaningful for version 5
    std::string title;
    std::string lmdb_file; // version 5 only
    std::string date;
    uint32_t num_oids = 0;
    uint64_t total_length = 0;
    uint32_t max_length = 0;
    std::vector<uint32_t> header_index;
    std::vector<uint32_t> sequence_index;
    size_t ambiguity_offsets_offset = 0; // nucleotide only
    size_t pin_length = 0;
};

std::string build_title(const std::vector<BlastDefLine>& deflines, const char* delimiter, bool all);
std::string format_seqid(const SeqId& id);

class BlastVolume {
public:

    struct RawChunk : public ::RawChunk {
        virtual bool empty() const override {
            return end_ <= begin_;
        }
        virtual DecodedPackage* decode(SequenceFile::Flags flags, const BitVector* filter, std::unordered_map<std::string, bool>* accs) const override;
        virtual OId begin() const noexcept override {
            return begin_;
        }
        virtual OId end() const noexcept override {
            return end_;
        }
        virtual size_t letters() const noexcept override {
            return letters_;
		}
        virtual size_t bytes() const noexcept override {
			return seq_data.size() + phr_data.size();
        }
		virtual ~RawChunk() override = default;
        std::vector<char> seq_data, phr_data;
		std::vector<uint32_t> seq_index, phr_index;
        OId begin_, end_;
        size_t letters_;
    };

    BlastVolume(const std::string& path, int idx, OId begin, OId end, bool load_index);
    const PinIndex& index() const { return index_; }
    std::vector<BlastDefLine> deflines(uint32_t oid, bool all, bool full_titles, bool taxids);
    std::vector<Letter> sequence(uint32_t oid);
    std::vector<char> raw_sequence(uint32_t count);
    std::vector<char> raw_deflines(uint32_t count);
	Loc length(uint32_t oid) const noexcept;
    size_t id_len(uint32_t oid) const noexcept;
    RawChunk* raw_chunk(size_t letters, SequenceFile::Flags flags);
	uint32_t seq_ptr() const noexcept { return seq_ptr_; }
    void rewind();
    int idx;
    OId begin = 0, end = 0;    

private:
    BlastVolume(PinIndex index, File&& phr, File&& psq, int idx, OId begin, OId end);
    static PinIndex ParsePinFile(File& mapping, bool index);
    
    PinIndex index_;
    File phr_mapping_;
    File psq_mapping_;
    uint32_t seq_ptr_ = 0, hdr_ptr_ = 0;
    
};

std::vector<BlastDefLine> decode_deflines(const char* header_data, size_t len, bool all, bool full_titles, bool taxids);
std::vector<Letter> decode_protein_sequence(const char* data, size_t len);