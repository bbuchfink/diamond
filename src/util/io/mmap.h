/****
Copyright © 2013-2025 Benjamin J. Buchfink

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
****/

#pragma once

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>
#include <filesystem>
#include <algorithm>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#else
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#ifdef __linux__
#ifndef O_TMPFILE
#define O_TMPFILE 020000000
#endif
#endif
#endif

#include <iostream>

class MMap {
public:
    struct Options {
        Options():
            initial_capacity(128ull * 1024),
			growth_factor(2.0)
        {}
        std::uint64_t initial_capacity = 128ull * 1024;
        double        growth_factor = 2.0;
        std::filesystem::path temp_dir_override;
    };

    explicit MMap(const Options& opts = Options()) {
        m_page = detect_page_size();
        m_growth = opts.growth_factor < 1.2 ? 1.2 : opts.growth_factor;
        create_temp_file(opts.temp_dir_override);
        map_resize(align_up(std::max<std::uint64_t>(opts.initial_capacity, m_page), m_page));
        m_size = 0;
    }

    MMap(const MMap&) = delete;
    MMap& operator=(const MMap&) = delete;

    MMap(MMap&& o) noexcept { move_from(std::move(o)); }
    MMap& operator=(MMap&& o) noexcept {
        if (this != &o) { close(); move_from(std::move(o)); }
        return *this;
    }

    ~MMap() { close(); }

    std::uint64_t write(const void* data, std::size_t bytes) {
        const std::uint64_t off = m_size;
        write_at(off, data, bytes);
        return off;
    }
    std::uint64_t write(std::string_view s) { return write(s.data(), s.size()); }

    void write_at(std::uint64_t offset, const void* data, std::size_t bytes) {
        if (bytes == 0) return;
        const std::uint64_t end = offset + bytes;
        ensure_capacity(end);
        std::memcpy(m_base + offset, data, bytes);
        if (end > m_size) m_size = end;
    }

    void reserve(std::uint64_t cap) {
        if (cap > m_capacity) map_resize(align_up(cap, m_page));
    }

    std::uint8_t* data() { return m_base; }
    const std::uint8_t* data() const { return m_base; }
    std::uint64_t size() const { return m_size; }
    std::uint64_t capacity() const { return m_capacity; }
    std::uint64_t page_size() const { return m_page; }

    void set_logical_size(std::uint64_t new_size) {
        if (new_size > m_capacity) throw std::runtime_error("logical size exceeds capacity");
        m_size = new_size;
    }

    void close() {
#ifdef _WIN32
        if (m_base) { UnmapViewOfFile(m_base); m_base = nullptr; }
        if (m_mapping) { CloseHandle(m_mapping);  m_mapping = nullptr; }
        if (m_file && m_file != INVALID_HANDLE_VALUE) { CloseHandle(m_file); m_file = INVALID_HANDLE_VALUE; }
#else
        if (m_base) { ::munmap(m_base, static_cast<size_t>(m_capacity)); m_base = nullptr; }
        if (m_fd >= 0) { ::close(m_fd); m_fd = -1; }
        if (!m_posix_path.empty()) {
            std::error_code ec;
            std::filesystem::remove(m_posix_path, ec);
            m_posix_path.clear();
        }
#endif
        m_capacity = 0; m_size = 0;
    }

private:
#ifdef _WIN32
    HANDLE m_file = INVALID_HANDLE_VALUE;
    HANDLE m_mapping = nullptr;
#else
    int m_fd = -1;
    std::filesystem::path m_posix_path;
#endif
    std::uint8_t* m_base = nullptr;
    std::uint64_t  m_size = 0;
    std::uint64_t  m_capacity = 0;
    std::uint64_t  m_page = 4096;
    double         m_growth = 2.0;

    static std::uint64_t align_up(std::uint64_t n, std::uint64_t a) {
        return (n + (a - 1)) / a * a;
    }

    std::uint64_t grow_strategy(std::uint64_t min_needed) const {
        std::uint64_t next = m_capacity > 0
            ? static_cast<std::uint64_t>(m_capacity * m_growth)
            : m_page * 4;
        if (next < min_needed) next = min_needed;
        return align_up(next, m_page);
    }

    void ensure_capacity(std::uint64_t required_end) {
        if (required_end <= m_capacity) return;
        map_resize(grow_strategy(required_end));
    }

    void map_resize(std::uint64_t new_capacity) {
#ifdef _WIN32
        if (m_base) { UnmapViewOfFile(m_base); m_base = nullptr; }
        if (m_mapping) { CloseHandle(m_mapping);  m_mapping = nullptr; }

        LARGE_INTEGER li; li.QuadPart = static_cast<LONGLONG>(new_capacity);
        if (!SetFilePointerEx(m_file, li, nullptr, FILE_BEGIN)) throw_win("SetFilePointerEx failed");
        if (!SetEndOfFile(m_file))                              throw_win("SetEndOfFile failed");

        DWORD hi = static_cast<DWORD>((new_capacity >> 32) & 0xFFFFFFFFull);
        DWORD lo = static_cast<DWORD>(new_capacity & 0xFFFFFFFFull);
        m_mapping = CreateFileMappingW(m_file, nullptr, PAGE_READWRITE, hi, lo, nullptr);
        if (!m_mapping) throw_win("CreateFileMappingW failed");

        void* p = MapViewOfFile(m_mapping, FILE_MAP_WRITE | FILE_MAP_READ, 0, 0, 0);
        if (!p) {
            DWORD err = GetLastError();
            CloseHandle(m_mapping); m_mapping = nullptr;
            throw_win("MapViewOfFile failed", err);
        }
        m_base = static_cast<std::uint8_t*>(p);
#else
        if (m_base) { ::munmap(m_base, static_cast<size_t>(m_capacity)); m_base = nullptr; }
        if (::ftruncate(m_fd, static_cast<off_t>(new_capacity)) != 0) throw_posix("ftruncate failed");

        void* p = ::mmap(nullptr, static_cast<size_t>(new_capacity),
            PROT_READ | PROT_WRITE, MAP_SHARED, m_fd, 0);
        if (p == MAP_FAILED) throw_posix("mmap failed");
        m_base = static_cast<std::uint8_t*>(p);
#endif
        m_capacity = new_capacity;
    }

    void create_temp_file(const std::filesystem::path& dir_override) {
#ifdef _WIN32
		wchar_t dir_buf[MAX_PATH];
        DWORD n = GetTempPathW(MAX_PATH, dir_buf);
        if (n == 0 || n > MAX_PATH) throw_win("GetTempPathW failed");
        GetCurrentDirectoryW(MAX_PATH, dir_buf);
        std::wstring dir = dir_buf;

        if (!dir_override.empty()) {
            dir = dir_override.wstring();
            if (!dir.empty() && dir.back() != L'\\' && dir.back() != L'/') dir.push_back(L'\\');
        }

        wchar_t name_buf[MAX_PATH];
        if (!GetTempFileNameW(dir.c_str(), L"diamond-tmp-", 0, name_buf))
            throw_win("GetTempFileNameW failed");
		std::wcout << name_buf << std::endl;

        m_file = CreateFileW(
            name_buf,
            GENERIC_READ | GENERIC_WRITE,
            FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
            nullptr,
            CREATE_ALWAYS,
            FILE_ATTRIBUTE_TEMPORARY | FILE_FLAG_DELETE_ON_CLOSE,
            nullptr
        );
        if (m_file == INVALID_HANDLE_VALUE) throw_win("CreateFileW failed");

#else
        std::filesystem::path base = dir_override; // dir_override.empty()
            //? std::filesystem::temp_directory_path()
            //: dir_override;

        bool used_otmpfile = false;
#ifdef __linux__
        int tmp_fd = ::open(base.c_str(),
            O_RDWR
#ifdef O_CLOEXEC
            | O_CLOEXEC
#endif
            | O_TMPFILE,
            S_IRUSR | S_IWUSR); // 0600
        if (tmp_fd >= 0) {
            m_fd = tmp_fd;
            used_otmpfile = true;
        }
#endif

        if (!used_otmpfile) {
            std::string tmpl = (base / "diamond-tmp-XXXXXX").string();
            std::vector<char> pathbuf(tmpl.begin(), tmpl.end());
            pathbuf.push_back('\0');

            m_fd = ::mkstemp(pathbuf.data());
            if (m_fd < 0) throw_posix("mkstemp failed");

            int flags = ::fcntl(m_fd, F_GETFD);
            if (flags >= 0) ::fcntl(m_fd, F_SETFD, flags | FD_CLOEXEC);

            m_posix_path = std::filesystem::path(pathbuf.data());
        }
#endif
    }

    std::uint64_t detect_page_size() {
#ifdef _WIN32
        SYSTEM_INFO si{}; GetSystemInfo(&si);
        return static_cast<std::uint64_t>(si.dwAllocationGranularity);
#else
        long p = ::sysconf(_SC_PAGESIZE);
        if (p <= 0) p = 4096;
        return static_cast<std::uint64_t>(p);
#endif
    }

#ifndef _WIN32
    [[noreturn]] static void throw_posix(const char* what) {
        int e = errno;
        std::string msg = std::string(what) + ": " + std::strerror(e);
        throw std::runtime_error(msg);
    }
#else
    [[noreturn]] static void throw_win(const char* what, DWORD code = GetLastError()) {
        LPVOID buf = nullptr;
        DWORD len = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
            FORMAT_MESSAGE_IGNORE_INSERTS,
            NULL, code, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
            (LPSTR)&buf, 0, NULL);
        std::string sys = (len && buf) ? std::string((LPSTR)buf, len) : std::string("Unknown error");
        if (buf) LocalFree(buf);
        std::string msg = std::string(what) + " (code " + std::to_string(code) + "): " + sys;
        throw std::runtime_error(msg);
    }
#endif

    void move_from(MMap&& o) noexcept {
        m_base = o.m_base;      o.m_base = nullptr;
        m_size = o.m_size;      o.m_size = 0;
        m_capacity = o.m_capacity; o.m_capacity = 0;
        m_page = o.m_page;
        m_growth = o.m_growth;
#ifdef _WIN32
        m_file = o.m_file;      o.m_file = INVALID_HANDLE_VALUE;
        m_mapping = o.m_mapping; o.m_mapping = nullptr;
#else
        m_fd = o.m_fd;          o.m_fd = -1;
        m_posix_path = std::move(o.m_posix_path);
        o.m_posix_path.clear();
#endif
    }
};