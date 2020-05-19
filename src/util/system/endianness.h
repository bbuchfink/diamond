#ifndef ENDIANNESS_H_
#define ENDIANNESS_H_

#include <stdint.h>

#ifdef WIN32

static inline uint64_t htole64(uint64_t host_64bits) {
	return host_64bits;
}

static inline uint32_t htole32(uint32_t host_32bits) {
	return host_32bits;
}

static inline uint32_t le32toh(uint32_t little_endian_32bits) {
	return little_endian_32bits;
}

static inline uint64_t le64toh(uint64_t little_endian_64bits) {
	return little_endian_64bits;
}

#else
#include <endian.h>
#endif

template<typename _t>
void to_host_endianness(_t &x) {}

template<>
static inline void to_host_endianness<uint64_t>(uint64_t& x)
{
	x = le64toh(x);
}

template<>
static inline void to_host_endianness<uint32_t>(uint32_t& x)
{
	x = le32toh(x);
}

template<typename _t>
_t to_little_endianness(_t x) {
	return x;
}

template<>
static inline uint64_t to_little_endianness<uint64_t>(uint64_t x)
{
	return htole64(x);
}

template<>
static inline uint32_t to_little_endianness<uint32_t>(uint32_t x)
{
	return htole32(x);
}

#endif