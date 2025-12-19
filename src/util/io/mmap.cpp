/****
Copyright  2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

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

#ifndef _WIN32
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "mmap.h"

MappedFile& MappedFile::operator=(MappedFile&& other) noexcept
{
    if (this == &other) {
        return *this;
    }

    unmap();

    data_ = other.data_;
    size_ = other.size_;
#if defined(_WIN32)
    hFile = other.hFile;
    hMap = other.hMap;
#endif

    other.data_ = nullptr;
    other.size_ = 0;
#if defined(_WIN32)
    other.hFile = INVALID_HANDLE_VALUE;
    other.hMap = nullptr;
#endif

    return *this;
}

bool MappedFile::map(const char* path)
{
#if defined(_WIN32)
    HANDLE h = CreateFileA(path, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (h == INVALID_HANDLE_VALUE) return false;
    LARGE_INTEGER sz;
    if (!GetFileSizeEx(h, &sz)) { CloseHandle(h); return false; }
    HANDLE hm = CreateFileMappingA(h, NULL, PAGE_READONLY, sz.HighPart, sz.LowPart, NULL);
    if (!hm) { CloseHandle(h); return false; }
    LPVOID addr = MapViewOfFile(hm, FILE_MAP_READ, 0, 0, 0);
    if (!addr) { CloseHandle(hm); CloseHandle(h); return false; }
    data_ = (const uint8_t*)addr;
    size_ = (size_t)sz.QuadPart;
    hFile = h; hMap = hm;
    return true;
#else
    int fd = ::open(path, O_RDONLY);
    if (fd < 0) return false;
    struct stat st;
    if (fstat(fd, &st) != 0) { ::close(fd); return false; }
    if (st.st_size <= 0) { ::close(fd); return false; }
    void* addr = mmap(nullptr, (size_t)st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
    if (addr == MAP_FAILED) { ::close(fd); return false; }
    data_ = (const uint8_t*)addr;
    size_ = (size_t)st.st_size;
    ::close(fd);
    return true;
#endif
}

void MappedFile::unmap()
{
#if defined(_WIN32)
    if (data_) UnmapViewOfFile((LPCVOID)data_);
    if (hMap) CloseHandle(hMap);
    if (hFile != INVALID_HANDLE_VALUE) CloseHandle(hFile);
#else
    if (data_ && size_) munmap((void*)data_, size_);
#endif
    data_ = nullptr; size_ = 0;
#if defined(_WIN32)
    hFile = INVALID_HANDLE_VALUE;
    hMap = nullptr;
#endif
}
