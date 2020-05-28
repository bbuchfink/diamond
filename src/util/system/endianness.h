#ifndef ENDIANNESS_H_
#define ENDIANNESS_H_

#include <stdint.h>
#include "../../lib/psnip/endian.h"

template<typename _t>
_t big_endian_byteswap(_t x) {
	return x;
}

template<>
inline uint64_t big_endian_byteswap<uint64_t>(uint64_t x)
{
	return psnip_endian_le64(x);
}

template<>
inline uint32_t big_endian_byteswap<uint32_t>(uint32_t x)
{
	return psnip_endian_le32(x);
}

#endif