/****
DIAMOND protein aligner
Copyright (C) 2019-2024 Max Planck Society for the Advancement of Science e.V.

Code developed by Klaus Reuter

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

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
#ifdef WIN32
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
    locked(false),
    file_name_(file_name)
{
#ifdef WIN32
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
#ifdef WIN32
    CloseHandle(hFile);
#else
    close(fd);
#endif
}

int FileStack::lock() {
    mtx_.lock();
#ifdef WIN32
    OVERLAPPED overlapvar;
    overlapvar.Offset = 0;
    overlapvar.OffsetHigh = 0;
    if(!LockFileEx(hFile, LOCKFILE_EXCLUSIVE_LOCK, 0, MAXDWORD, MAXDWORD, &overlapvar))
        throw std::runtime_error("could not put lock on file " + file_name_);
    locked = true;
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
        } else {
            locked = true;
        }
    } else {
        throw(std::runtime_error("could not put lock on non-open file " + file_name_));
    }
    return fcntl_status;
#endif
}

int FileStack::unlock() {
#ifdef WIN32
    OVERLAPPED overlapvar;
    overlapvar.Offset = 0;
    overlapvar.OffsetHigh = 0;
    if (!UnlockFileEx(hFile, 0, MAXDWORD, MAXDWORD, &overlapvar))
        throw(std::runtime_error("could not unlock file " + file_name_));
    locked = false;
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
        } else {
            locked = false;
        }
    }
    mtx_.unlock();
    return fcntl_status;
#endif
}

int64_t FileStack::seek(int64_t offset, int mode) {
#ifdef WIN32
    const DWORD off = SetFilePointer(hFile, offset, NULL, mode == SEEK_END ? FILE_END : FILE_BEGIN);
    if (off == INVALID_SET_FILE_POINTER)
        throw std::runtime_error("Failed to seek");
    return off;
#else
    return lseek(fd, offset, mode);
#endif
}

int64_t FileStack::read(char* buf, int64_t size) {
#ifdef WIN32
    DWORD n;
    if (!ReadFile(hFile, buf, size, &n, NULL))
        throw std::runtime_error("Error reading file " + file_name_);
    return n;
#else
    return ::read(fd, buf, size);
#endif
}

int64_t FileStack::write(const char* buf, int64_t size) {
#ifdef WIN32
    DWORD n;
    if (!WriteFile(hFile, buf, size, &n, NULL))
        throw std::runtime_error("Error writing file " + file_name_);
    return n;
#else
    return ::write(fd, buf, size);
#endif
}

int FileStack::truncate(int64_t size) {
#ifdef WIN32
    seek(size, SEEK_SET);
    if (!SetEndOfFile(hFile))
        throw std::runtime_error("Error calling SetEndOfFile");
    return 0;
#else
    return ftruncate(fd, size);
#endif
}

int FileStack::pop_non_locked(string & buf, const bool keep_flag, size_t & size_after_pop) {
    DBG("");
    buf.clear();
    int stat = 0;
    const off_t size = seek(0, SEEK_END);
    if (size > 0) {
        off_t jmp = size - max_line_length;
        if (jmp < 0) jmp = 0;
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

int FileStack::pop_non_locked(string & buf) {
    DBG("");
    size_t size_after_pop = numeric_limits<size_t>::max();
    return pop_non_locked(buf, false, size_after_pop);
}

int FileStack::pop(string & buf, const bool keep_flag, size_t & size_after_pop) {
    DBG("");
    bool locked_internally = false;
    if (! locked) {
        lock();
        locked_internally = true;
    }
    int val = pop_non_locked(buf, keep_flag, size_after_pop);
    if (locked_internally) {
        unlock();
    }
    return val;
}

int FileStack::pop(string & buf) {
    DBG("");
    size_t size_after_pop = numeric_limits<size_t>::max();
    return pop(buf, false, size_after_pop);
}

int FileStack::top(string & buf) {
    DBG("");
    size_t size_after_pop = numeric_limits<size_t>::max();
    return pop(buf, true, size_after_pop);
}

int FileStack::pop(string & buf, size_t & size_after_pop) {
    DBG("");
    return pop(buf, false, size_after_pop);
}

int FileStack::pop(int & i) {
    DBG("");
    string buf;
    size_t size_after_pop = numeric_limits<size_t>::max();
    const int get_status = pop(buf, false, size_after_pop);
    if (get_status > 0) {
        return i = stoi(buf);
    } else {
        return i = -1;
    }
}

int FileStack::top(int & i) {
    DBG("");
    string buf;
    size_t size_after_pop = numeric_limits<size_t>::max();
    const int get_status = pop(buf, true, size_after_pop);
    if (get_status > 0) {
        return i = stoi(buf);
    } else {
        return i = -1;
    }
}

void FileStack::remove(const string & line) {
    DBG("");
    bool locked_internally = false;
    if (! locked) {
        lock();
        locked_internally = true;
    }
    const off_t size = seek(0, SEEK_END);
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

    if (locked_internally) {
        unlock();
    }
}

int64_t FileStack::push_non_locked(const string & buf) {
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
    bool locked_internally = false;
    if (! locked) {
        lock();
        locked_internally = true;
    }
    int64_t n = push_non_locked(buf);
    if (size_after_push != numeric_limits<size_t>::max()) {
        size_after_push = size();
    }
    if (locked_internally) {
        unlock();
    }
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

int64_t FileStack::push(int i) {
    DBG("");
    string buf = to_string(i);
    return push(buf);
}

size_t FileStack::size() {
    DBG("");
    size_t c = 0;
    const size_t chunk_size = default_max_line_length;
    char * raw = new char[chunk_size * sizeof(char)];

    bool locked_internally = false;
    if (! locked) {
        lock();
        locked_internally = true;
    }

    size_t n_bytes, i;
    seek(0, SEEK_SET);
    while ((n_bytes = read(raw, chunk_size)) > 0) {
        for (i=0; i<n_bytes; i++) {
            if (raw[i] == '\n') {
                c++;
            }
        }
    }

    if (locked_internally) {
        unlock();
    }

    delete [] raw;
    return c;
}

int FileStack::clear() {
    DBG("");

    bool locked_internally = false;
    if (! locked) {
        lock();
        locked_internally = true;
    }

    seek(0, SEEK_SET);
    int stat = truncate(0);

    if (locked_internally) {
        unlock();
    }

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