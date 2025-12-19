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

#include <assert.h>
#include "volume.h"
#include "ber.h"

using ByteView = MappedFile::View;
using std::runtime_error;

PinIndex Volume::ParsePinFile(const MappedFile& mapping)
{
    const ByteView data = mapping.view();
    std::size_t offset = 0;

    PinIndex index;
    index.version = ReadBE32(data, offset);
    if (index.version != 4 && index.version != 5) {
        throw runtime_error("Unsupported database format version: " + std::to_string(index.version));
    }

    const std::uint32_t seq_type_flag = ReadBE32(data, offset);
    index.is_protein = (seq_type_flag == 1);

    if (index.version == 5) {
        index.volume_number = ReadBE32(data, offset);
    }

    index.title = ReadPascalString(data, offset);

    if (index.version == 5) {
        index.lmdb_file = ReadPascalString(data, offset);
    }

    index.date = ReadPascalString(data, offset);

    index.num_oids = ReadBE32(data, offset);
    index.total_length = ReadLE64(data, offset);
    index.max_length = ReadBE32(data, offset);
    index.pin_length = data.size();

    const std::size_t count = static_cast<std::size_t>(index.num_oids) + 1;
    const auto offsets_size = count * sizeof(std::uint32_t);
    auto reserve_offset_array = [&](std::size_t& target) {
        if (offset + offsets_size > data.size()) {
            throw runtime_error("Offset array exceeds PIN file size");
        }
        target = offset;
        offset += offsets_size;
        };

    reserve_offset_array(index.header_offsets_offset);
    reserve_offset_array(index.sequence_offsets_offset);

    if (!index.is_protein) {
        reserve_offset_array(index.ambiguity_offsets_offset);
    }

    return index;
}