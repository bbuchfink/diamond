/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef PACKED_LOC_H_
#define PACKED_LOC_H_

struct packed_uint40_t
{
	uint8_t		high;
	uint32_t	low;
	packed_uint40_t():
		high (),
		low ()
	{ }
	packed_uint40_t(uint64_t v):
		high (v>>32),
		low (v&0xffffffffu)
	{ }
	operator const uint64_t() const
	{ return (uint64_t(high) << 32) | low; }
	bool operator<(const packed_uint40_t &rhs) const
	{ return high < rhs.high || (high == rhs.high && low < rhs.low); }
	friend uint64_t operator-(const packed_uint40_t &x, const packed_uint40_t &y)
	{ return (const uint64_t)(x) - (const uint64_t)(y); }
} __attribute__((packed));

template<class _loc>
struct packed_sequence_location
{
	typedef _loc type;
};

template<>
struct packed_sequence_location<uint64_t>
{
	typedef packed_uint40_t type;
};

#endif /* PACKED_LOC_H_ */
