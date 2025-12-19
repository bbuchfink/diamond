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
#include "util/system.h"

using std::vector;
using ByteView = MappedFile::View;
using std::out_of_range;
using std::runtime_error;

vector<Letter> DecodeProteinSequence(const ByteView& data, std::uint32_t start, std::uint32_t end)
{
    if (start > end || end > data.size()) {
        throw out_of_range("Sequence offsets exceed sequence file length");
    }

    vector<Letter> decoded;
    decoded.reserve(end - start);
    for (std::uint32_t pos = start; pos < end; ++pos) {
        const uint8_t aa = data[pos];
        if (aa == '\0') {
            break;
        }
        if (aa >= sizeof(NCBI_TO_STD) / sizeof(NCBI_TO_STD[0]))
            hard_fail("Invalid amino acid code in sequence data");
        decoded.push_back(NCBI_TO_STD[aa]);
    }
    return decoded;
}

vector<Letter> Volume::sequence(std::uint32_t oid) const
{
    if (oid >= index_.num_oids) {
        throw out_of_range("OID exceeds number of sequences in volume");
    }

    const ByteView pin_view = pin_mapping_.view();
    const ByteView psq_view = psq_mapping_.view();
    const std::uint32_t start = read_offset(pin_view, index_.sequence_offsets_offset, oid);
    const std::uint32_t end = read_offset(pin_view, index_.sequence_offsets_offset, oid + 1);

    if (!index_.is_protein) {
        throw runtime_error("Nucleotide sequence decoding is not supported yet");
    }

    return DecodeProteinSequence(psq_view, start, end);
}

Loc Volume::length(uint32_t oid) const {
    const ByteView pin_view = pin_mapping_.view();
    return read_offset(pin_view, index_.sequence_offsets_offset, oid + 1) - read_offset(pin_view, index_.sequence_offsets_offset, oid) - 1;
}