#ifndef ENDIANNESS_H_
#define ENDIANNESS_H_

#include <stdint.h>
#include "../../lib/psnip/endian.h"

template<typename _t>
void to_host_endianness(_t &x) {}

template<>
inline void to_host_endianness<uint64_t>(uint64_t& x)
{
	x = psnip_endian_le64(x);
}

template<>
inline void to_host_endianness<uint32_t>(uint32_t& x)
{
	x = psnip_endian_le32(x);
}

template<typename _t>
_t to_little_endianness(_t x) {
	return x;
}

template<>
inline uint64_t to_little_endianness<uint64_t>(uint64_t x)
{
	return psnip_endian_le64(x);
}

template<>
inline uint32_t to_little_endianness<uint32_t>(uint32_t x)
{
	return psnip_endian_le32(x);
}

#endif