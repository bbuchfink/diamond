#pragma once

namespace mio {

enum class access_mode
{
    read,
    write
};

template<access_mode AccessMode, typename ByteT>
struct basic_mmap;

template<typename ByteT>
using basic_mmap_source = basic_mmap<access_mode::read, ByteT>;

using mmap_source = basic_mmap_source<char>;

}