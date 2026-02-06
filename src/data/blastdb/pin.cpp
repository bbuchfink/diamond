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

#include <set>
#include <assert.h>
#include "volume.h"
#include "ber.h"

using std::set;
using std::vector;
using std::runtime_error;
using std::string;
using std::unordered_map;

PinIndex BlastVolume::ParsePinFile(File& mapping, bool load_index)
{
    size_t offset = 0;

    PinIndex index;    
    index.version = ReadBE32(mapping);
    if (index.version != 4 && index.version != 5) {
        throw runtime_error("Unsupported database format version: " + std::to_string(index.version));
    }

    const uint32_t seq_type_flag = ReadBE32(mapping);
    index.is_protein = (seq_type_flag == 1);

    if (index.version == 5)
        index.volume_number = ReadBE32(mapping);
    index.title = ReadPascalString(mapping);

    if (index.version == 5)
        index.lmdb_file = ReadPascalString(mapping);

    index.date = ReadPascalString(mapping);
    index.num_oids = ReadBE32(mapping);    
    index.total_length = ReadLE64(mapping);    
    index.max_length = ReadBE32(mapping);
    //index.pin_length = data.size();
    if (!load_index)
        return index;

    const size_t count = static_cast<std::size_t>(index.num_oids) + 1;
	index.header_index.reserve(count);
	index.sequence_index.reserve(count);
    for (size_t i = 0; i < count; ++i)
        index.header_index.push_back(ReadBE32(mapping));

    for (size_t i = 0; i < count; ++i)
        index.sequence_index.push_back(ReadBE32(mapping));

    if (!index.is_protein) {
        //reserve_offset_array(index.ambiguity_offsets_offset);
    }

    return index;
}

BlastVolume::BlastVolume(const string& path, int idx, OId begin, OId end, bool load_index) :
    idx(idx),
    begin(begin),
    end(end),
    phr_mapping_(path + ".phr", "rb"),
    psq_mapping_(path + ".psq", "rb")    
{
	File pin(path + ".pin", "rb");
    index_ = BlastVolume::ParsePinFile(pin, load_index);
    pin.close();        
}

BlastVolume::RawChunk* BlastVolume::raw_chunk(size_t letters, SequenceFile::Flags flags) {
    uint32_t begin = hdr_ptr_;
    if (!bool(flags & SequenceFile::Flags::SEQS)) {
        if(seq_ptr_ != 0)
			throw runtime_error("Volume::raw_chunk");
    } else if (!bool((flags & SequenceFile::Flags::TITLES) | (flags & SequenceFile::Flags::TAXON_MAPPING))) {
        if (hdr_ptr_ != 0)
            throw runtime_error("Volume::raw_chunk");
		begin = seq_ptr_;
    }
    else if (hdr_ptr_ != seq_ptr_)
        throw runtime_error("Cannot read raw chunk: last accessed header and sequence OIDs do not match");
    uint32_t end = begin;
    size_t l = 0;
    while (end < index_.num_oids && l < letters) {
        l += length(end);
        ++end;
    }
    RawChunk* chunk = new RawChunk();
    chunk->letters_ = l;
	chunk->begin_ = begin + this->begin;
    chunk->end_ = end + this->begin;
    const uint32_t n = end - begin;
    if (n == 0)
        return chunk;
    if (bool((flags & SequenceFile::Flags::TITLES) | (flags & SequenceFile::Flags::TAXON_MAPPING))) {
        chunk->phr_index.assign(index_.header_index.begin() + hdr_ptr_, index_.header_index.begin() + hdr_ptr_ + n + 1);
        chunk->phr_data = raw_deflines(n);
    }
    if (bool(flags & SequenceFile::Flags::SEQS)) {
        chunk->seq_index.assign(index_.sequence_index.begin() + seq_ptr_, index_.sequence_index.begin() + seq_ptr_ + n + 1);
        chunk->seq_data = raw_sequence(n);
    }
    return chunk;
}

static bool acc_filter(const vector<BlastDefLine>& deflines, unordered_map<std::string, bool>& accs) {
    for (const auto& d : deflines) {
        for (const auto& s : d.seqids) {
            auto it = accs.find(s.value);
            if (it == accs.end() && (s.version || s.chain))
                it = accs.find(format_seqid(s));
            if (it != accs.end()) {
                it->second = true;
                return true;
            }
        }
    }
    return false;
}

DecodedPackage* BlastVolume::RawChunk::decode(SequenceFile::Flags flags, const BitVector* filter, unordered_map<std::string, bool>* accs) const {
    assert(filter == nullptr || accs == nullptr);
    DecodedPackage* pkg = new DecodedPackage();
    pkg->no = no;
	const size_t n = end_ - begin_;
    const char* seq_ptr = seq_data.data(), *phr_ptr = phr_data.data();
    const bool titles = bool(flags & SequenceFile::Flags::TITLES), seqs = bool(flags & SequenceFile::Flags::SEQS), taxids = bool(flags & SequenceFile::Flags::TAXON_MAPPING),
        full_titles = bool(flags & SequenceFile::Flags::FULL_TITLES), all_seqids = bool(flags & SequenceFile::Flags::ALL_SEQIDS);
    pkg->oids.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        const OId oid = begin_ + i;
        bool f = !filter || filter->get(oid);
        
        if (titles || taxids || accs) {
            const size_t lhdr = phr_index[i + 1] - phr_index[i];
            if (f || accs) {
                const vector<BlastDefLine> deflines = decode_deflines(phr_ptr, lhdr, all_seqids, full_titles, taxids);
                if (accs)
                    f = acc_filter(deflines, *accs);
                if (f && titles) {
                    const string title = build_title(deflines, "\1", true);
                    pkg->ids.push_back(title.begin(), title.end());
                }
                if (f && taxids) {
                    set<TaxId> s;
                    for (auto j = deflines.begin(); j != deflines.end(); ++j)
                        if (j->taxid)
                            s.insert(j->taxid.value());
                    for (TaxId t : s)
                        pkg->taxids.emplace_back(oid, t);
                }
            }
            phr_ptr += lhdr;
        }
           
        if (seqs) {
            const size_t lseq = seq_index[i + 1] - seq_index[i];
            if (f) {
                const vector<Letter> seq = decode_protein_sequence(seq_ptr, lseq);
                pkg->seqs.push_back(seq.begin(), seq.end());
            }
            seq_ptr += lseq;
        }

        if (f)
            pkg->oids.push_back(oid);
    }
    return pkg;
}

void BlastVolume::rewind() {
    hdr_ptr_ = 0;
    seq_ptr_ = 0;
	phr_mapping_.seek(0, SEEK_SET);
	psq_mapping_.seek(0, SEEK_SET);
}