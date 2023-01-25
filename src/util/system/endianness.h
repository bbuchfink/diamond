#ifndef ENDIANNESS_H_
#define ENDIANNESS_H_

#include <stdint.h>
#include "../../lib/psnip/endian.h"

static inline bool is_little_endian() {
	static int32_t test = 1;
	return *reinterpret_cast<int8_t*>(&test) == 1;
}

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

template<>
inline uint16_t big_endian_byteswap<uint16_t>(uint16_t x)
{
	return psnip_endian_le16(x);
}

template<>
inline int64_t big_endian_byteswap<int64_t>(int64_t x)
{
	return psnip_endian_le64(x);
}

template<>
inline int32_t big_endian_byteswap<int32_t>(int32_t x)
{
	return psnip_endian_le32(x);
}

template<>
inline int16_t big_endian_byteswap<int16_t>(int16_t x)
{
	return psnip_endian_le16(x);
}

#endif