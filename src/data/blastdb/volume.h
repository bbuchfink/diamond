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
	Loc length(uint32_t oid);
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