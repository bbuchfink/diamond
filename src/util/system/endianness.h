#ifndef ENDIANNESS_H_
#define ENDIANNESS_H_

#ifdef WIN32

#include <stdint.h>

uint64_t htole64(uint64_t host_64bits) {
	return host_64bits;
}

#else
#include <endian.h>
#endif

#endif