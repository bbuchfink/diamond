/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Code developed by Klaus Reuter

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#include <string>
#include <iostream>
#include <chrono>
#include <thread>
#include <stdexcept>
#include <limits>
#include <vector>
#include <algorithm>
#include <string.h>
#include "multiprocessing.h"
#ifdef _WIN32
#else
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

// #define DEBUG
#undef DEBUG
#include "filestack.h"

using std::runtime_error;
using std::this_thread::sleep_for;
using std::numeric_limits;
using std::to_string;
using std::string;
using std::vector;

const string default_file_name = "default_stack.idx";
const int default_max_line_length = 4096;


FileStack::FileStack() : FileStack::FileStack(default_file_name) {
    std::cerr << "FileStack: Using default file name " << default_file_name << std::endl;
}

FileStack::FileStack(const string & file_name) : FileStack::FileStack(file_name, default_max_line_length) {
}

FileStack::FileStack(const string & file_name, int maximum_line_length):
    file_name_(file_name)
{
#ifdef _WIN32
    hFile = CreateFile(TEXT(file_name.c_str()), FILE_GENERIC_WRITE | FILE_GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, OPEN_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE)
        throw std::runtime_error("Error opening file " + file_name_);
#else
    DBG("");
    fd = open(file_name.c_str(), O_RDWR | O_CREAT, 00664);
    if (fd == -1) {
        throw(std::runtime_error("could not open file " + file_name_));
    }
#endif
    set_max_line_length(maximum_line_length);
}

FileStack::~FileStack() {
    DBG("");
#ifdef _WIN32
    CloseHandle(hFile);
#else
    close(fd);
#endif
}

int FileStack::lock() {
    mtx_.lock();
#ifdef _WIN32
    OVERLAPPED overlapvar;
    overlapvar.Offset = 0;
    overlapvar.OffsetHigh = 0;
    if(!LockFileEx(hFile, LOCKFILE_EXCLUSIVE_LOCK, 0, MAXDWORD, MAXDWORD, &overlapvar))
        throw runtime_error("could not put lock on file " + file_name_);
    return 0;
#else
    DBG("");
    int fcntl_status = -1;
    if (fd >= 0) {
        memset(&lck, 0, sizeof(lck));
        lck.l_type = F_WRLCK; // exclusive lock
        lck.l_whence = SEEK_SET;
        lck.l_start = 0;
        lck.l_len = 0;
        fcntl_status = fcntl(fd, F_SETLKW, &lck);
        if (fcntl_status == -1) {
            throw(std::runtime_error("could not put lock on file " + file_name_));
        }
    } else {
        throw(std::runtime_error("could not put lock on non-open file " + file_name_));
    }
    return fcntl_status;
#endif
}

int FileStack::unlock() {
#ifdef _WIN32
    OVERLAPPED overlapvar;
    overlapvar.Offset = 0;
    overlapvar.OffsetHigh = 0;
    if (!UnlockFileEx(hFile, 0, MAXDWORD, MAXDWORD, &overlapvar))
        throw(std::runtime_error("could not unlock file " + file_name_));
    mtx_.unlock();
    return 0;
#else
    DBG("");
    int fcntl_status = -1;
    if (fd >= 0) {
        lck.l_type = F_UNLCK;
        fcntl_status = fcntl(fd, F_SETLKW, &lck);
        if (fcntl_status == -1) {
            throw(std::runtime_error("could not unlock file " + file_name_));
        }
    }
    mtx_.unlock();
    return fcntl_status;
#endif
}

int64_t FileStack::seek(int64_t offset, int mode) {
#ifdef _WIN32
    LARGE_INTEGER ptr, o;
    o.QuadPart = offset;
    const BOOL ret = SetFilePointerEx(hFile, o, &ptr, mode == SEEK_END ? FILE_END : FILE_BEGIN);
    if (!ret)
        throw std::runtime_error("Failed to seek");
    return ptr.QuadPart;
#else
    return lseek(fd, offset, mode);
#endif
}

size_t FileStack::read(char* buf, size_t size) {
#ifdef _WIN32
    DWORD n;
    if (!ReadFile(hFile, buf, (DWORD)size, &n, NULL))
        throw runtime_error("Error reading file " + file_name_);
    return n;
#else
    ssize_t i = ::read(fd, buf, size);
    if (i < 0)
        throw runtime_error("Error reading file " + file_name_);
    return (size_t)i;
#endif
}

int64_t FileStack::write(const char* buf, size_t size) {
#ifdef _WIN32
    DWORD n;
    if (!WriteFile(hFile, buf, (DWORD)size, &n, NULL))
        throw std::runtime_error("Error writing file " + file_name_);
    return n;
#else
    return ::write(fd, buf, size);
#endif
}

int FileStack::truncate(size_t size) {
#ifdef _WIN32
    seek(size, SEEK_SET);
    if (!SetEndOfFile(hFile))
        throw std::runtime_error("Error calling SetEndOfFile");
    return 0;
#else
    return ftruncate(fd, size);
#endif
}

int64_t FileStack::pop_non_locked(string & buf, const bool keep_flag, size_t & size_after_pop) {
    DBG("");
    buf.clear();
    int stat = 0;
    const int64_t size = seek(0, SEEK_END);
    if (size > 0) {
        int64_t jmp = size >= max_line_length ? size - max_line_length : 0;
        seek(jmp, SEEK_SET);

        char * raw = new char[max_line_length * sizeof(char)];
        const int64_t n_read = read(raw, max_line_length);
        string chunk;
        chunk.assign(raw, n_read);
        delete [] raw;

        size_t begin = 0, end = 0;
        const char key = '\n';
        size_t found = chunk.rfind(key);
        if (found != string::npos) {
            end = found;
        }
        if (end > 0) {
            found = chunk.rfind(key, end-1);
            if (found != string::npos) {
                begin = found + 1;
            }
        }
        const size_t line_size = end - begin + 1;
        if (line_size > 0) {
            buf.assign(chunk, begin, line_size - (chunk[begin + line_size - 1] == '\n' ? 1 : 0));
            if (! keep_flag) {
                stat = truncate(size - line_size);
            }
        }
    }
    if (size_after_pop != numeric_limits<size_t>::max()) {
        size_after_pop = this->size();
    }
    if (stat == -1) {
        return stat;
    } else {
        DBG(buf);
        return buf.size();
    }
}

int64_t FileStack::pop_exclusive(string & buf) {
    DBG("");
    size_t size_after_pop = numeric_limits<size_t>::max();
    return pop_non_locked(buf, false, size_after_pop);
}

int FileStack::pop(string & buf, const bool keep_flag, size_t & size_after_pop) {
    DBG("");
    lock();
    int val = (int)pop_non_locked(buf, keep_flag, size_after_pop);
    unlock();
    return val;
}

int FileStack::pop(string& buf) {
    DBG("");
    size_t size_after_pop = numeric_limits<size_t>::max();
    return pop(buf, false, size_after_pop);
}

int FileStack::top(string & buf) {
    DBG("");
    size_t size_after_pop = numeric_limits<size_t>::max();
    return pop(buf, true, size_after_pop);
}

int FileStack::pop(string& buf, size_t& size_after_pop) {
    DBG("");
    return pop(buf, false, size_after_pop);
}

int64_t FileStack::pop_non_locked(int64_t & i) {
    DBG("");
    string buf;
    size_t size_after_pop = numeric_limits<size_t>::max();
    const int get_status = pop_non_locked(buf, false, size_after_pop);
    if (get_status > 0) {
        return i = std::stoll(buf);
    } else {
        return i = -1;
    }
}

int64_t FileStack::pop(int64_t& i) {
    lock();
    const int64_t r = pop_non_locked(i);
    unlock();
    return r;
}

int64_t FileStack::top(int64_t& i) {
    DBG("");
    string buf;
    size_t size_after_pop = numeric_limits<size_t>::max();
    const int get_status = pop(buf, true, size_after_pop);
    if (get_status > 0) {
        return i = std::stoll(buf);
    } else {
        return i = -1;
    }
}

void FileStack::remove(const string & line) {
    DBG("");
    lock();
    const int64_t size = seek(0, SEEK_END);
    seek(0, SEEK_SET);

    char * raw = new char[size * sizeof(char)];
    const int64_t n_read = read(raw, size);
    string buf;
    buf.assign(raw, n_read);
    delete [] raw;

    vector<string> tokens = split(buf, '\n');
    buf.clear();

    tokens.erase(std::remove(tokens.begin(), tokens.end(), line), tokens.end());

    seek(0, SEEK_SET);
    int stat = truncate(0);
    for (auto it = tokens.begin(); it != tokens.end(); ++it) {
        buf = *it + '\n';
        size_t n = write(buf.c_str(), buf.size());
    }

    unlock();
}

int64_t FileStack::push_exclusive(const string & buf) {
    DBG("");
    static const string nl("\n");
    seek(0, SEEK_END);
    size_t n = write(buf.c_str(), buf.size());
    if (buf.back() != nl.back()) {
        n += write(nl.c_str(), nl.size());
    }
    return n;
}

int64_t FileStack::push(const string & buf, size_t & size_after_push) {
    DBG("");
    lock();
    int64_t n = push_exclusive(buf);
    if (size_after_push != numeric_limits<size_t>::max()) {
        size_after_push = size();
    }
    unlock();
    return n;
}

std::string FileStack::file_name() const {
    return file_name_;
}

int64_t FileStack::push(const string & buf) {
    DBG("");
    size_t size_after_push = numeric_limits<size_t>::max();
    return push(buf, size_after_push);
}

int64_t FileStack::push_non_locked(int64_t i) {
    DBG("");
    string buf = to_string(i);
    return push_exclusive(buf);
}

int64_t FileStack::push(int64_t i) {
    lock();
    const int64_t r = push_non_locked(i);
    unlock();
    return r;
}

size_t FileStack::size() {
    DBG("");
    size_t c = 0;
    const size_t chunk_size = default_max_line_length;
    char * raw = new char[chunk_size * sizeof(char)];

    lock();

    size_t n_bytes, i;
    seek(0, SEEK_SET);
    while ((n_bytes = read(raw, chunk_size)) > 0) {
        for (i=0; i<n_bytes; i++) {
            if (raw[i] == '\n') {
                c++;
            }
        }
    }

    unlock();

    delete [] raw;
    return c;
}

int FileStack::clear() {
    DBG("");
    lock();
    seek(0, SEEK_SET);
    int stat = truncate(0);
    unlock();
    return stat;
}

bool FileStack::poll_query(const string & query, const double sleep_s, const size_t max_iter) {
    DBG("");
    string buf;
    const std::chrono::duration<double> sleep_time(sleep_s);
    for (size_t i=0; i < max_iter; ++i) {
        top(buf);
        if (buf.find(query) != string::npos) {
            DBG(string("") + " success : poll_iteration=" + to_string(i) + ", query=" + "\"" + query + "\"");
            return true;
        } else {
            DBG(string("") + " ongoing : poll_iteration=" + to_string(i) + ", query=" + "\"" + query + "\"");
        }
        if (buf.find("STOP") != string::npos) {
            throw(runtime_error("STOP on FileStack " + file_name_));
        }
        sleep_for(sleep_time);
    }
    throw(runtime_error("Could not discover keyword " + query + " on FileStack " + file_name_
                      + " within " + to_string(double(max_iter) * sleep_s) + " seconds."));
    return false;  // TODO : finally decide on the semantics
};

bool FileStack::poll_size(const size_t size, const double sleep_s, const size_t max_iter) {
    DBG("");
    string buf;
    const std::chrono::duration<double> sleep_time(sleep_s);
    for (size_t i=0; i < max_iter; ++i) {
        if (this->size() == size) {
            return true;
        }
        sleep_for(sleep_time);
    }
    throw(runtime_error("Could not detect size " + to_string(size) + " of FileStack " + file_name_
                      + " within " + to_string(double(max_iter) * sleep_s) + " seconds."));
    return false;  // TODO : finally decide on the semantics
};



int FileStack::set_max_line_length(int n) {
    DBG("");
    const int minimum_line_length = 8;
    if (n < minimum_line_length) {
        n = minimum_line_length;
    }
    max_line_length = n;
    return max_line_length;
}

int FileStack::get_max_line_length() {
    DBG("");
    return max_line_length;
}


int64_t FileStack::fetch_add(int64_t n) {
    string s;
    lock();
    int64_t i;
    pop_non_locked(i);
    if (i == -1)
        i = 0;
    push_non_locked(i + n);
    // message_stream << "Atomic = " << i << std::endl;
    unlock();
    return i;
}