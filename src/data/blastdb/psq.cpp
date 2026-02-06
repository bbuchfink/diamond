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

#include "volume.h"

using std::vector;
using std::out_of_range;
using std::runtime_error;

vector<Letter> decode_protein_sequence(const char* data, size_t len)
{
    vector<Letter> decoded;
    decoded.reserve(len);
    for (size_t i = 0; i < len; ++i) {
        const uint8_t aa = data[i];
        if (aa == '\0') {
            if (i == 0)
                continue;
            else if (i == len - 1)
                break;
            else
				throw runtime_error("Unexpected null terminator in sequence data");
        }
        if (aa >= sizeof(NCBI_TO_STD) / sizeof(NCBI_TO_STD[0]))
            throw runtime_error("Invalid amino acid code in sequence data");
        decoded.push_back(NCBI_TO_STD[aa]);
    }
    return decoded;
}

vector<Letter> BlastVolume::sequence(uint32_t oid)
{
    if (oid >= index_.num_oids) {
        throw out_of_range("OID exceeds number of sequences in volume");
    }

    const uint32_t start = index_.sequence_index[oid];
    const uint32_t end = index_.sequence_index[oid + 1];

    if (!index_.is_protein) {
        throw runtime_error("Nucleotide sequence decoding is not supported yet");
    }
    if (oid != seq_ptr_)
        psq_mapping_.seek(start, SEEK_SET);
    seq_ptr_ = oid + 1;
    return decode_protein_sequence(psq_mapping_.read(end - start), end - start);
}

vector<char> BlastVolume::raw_sequence(uint32_t count) {
    const size_t n = index_.sequence_index[seq_ptr_ + count] - index_.sequence_index[seq_ptr_];
    vector<char> v(n);
    psq_mapping_.read(v.data(), n);
    seq_ptr_ += count;
    return v;
}

Loc BlastVolume::length(uint32_t oid) const noexcept {
    return index_.sequence_index[oid + 1] - index_.sequence_index[oid] - 1;
}