#include <string>
#include <iostream>
#include <chrono>
#include <thread>
#include <stdexcept>
#include <limits>
#include <vector>
#include <algorithm>

#include <cstdio>
#include <cstring>
#ifndef WIN32
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "multiprocessing.h"
// #define DEBUG
#undef DEBUG
#include "filestack.h"

using namespace std;


const string default_file_name = "default_stack.idx";
const int default_max_line_length = 4096;


FileStack::FileStack() : FileStack::FileStack(default_file_name) {
    std::cerr << "FileStack: Using default file name " << default_file_name << std::endl;
}

FileStack::FileStack(const string & file_name) : FileStack::FileStack(file_name, default_max_line_length) {
}

FileStack::FileStack(const string & file_name, int maximum_line_length) {
#ifndef WIN32
    DBG("");
    fd = open(file_name.c_str(), O_RDWR | O_CREAT, 00664);
    if (fd == -1) {
        throw(std::runtime_error("could not open file " + file_name));
    }
    this->locked = false;
    this->file_name = file_name;
    set_max_line_length(maximum_line_length);
#endif
}

FileStack::~FileStack() {
    DBG("");
#ifndef WIN32
    close(fd);
#endif
}



int FileStack::lock() {
#ifndef WIN32
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
            throw(std::runtime_error("could not put lock on file " + file_name));
        } else {
            locked = true;
        }
    } else {
        throw(std::runtime_error("could not put lock on non-open file " + file_name));
    }
    return fcntl_status;
#else
    return 0;
#endif
}

int FileStack::unlock() {
#ifndef WIN32
    DBG("");
    int fcntl_status = -1;
    if (fd >= 0) {
        lck.l_type = F_UNLCK;
        fcntl_status = fcntl(fd, F_SETLKW, &lck);
        if (fcntl_status == -1) {
            throw(std::runtime_error("could not unlock file " + file_name));
        } else {
            locked = false;
        }
    }
    return fcntl_status;
#else
    return 0;
#endif
}

int FileStack::pop_non_locked(string & buf, const bool keep_flag, size_t & size_after_pop) {
#ifndef WIN32
    DBG("");
    buf.clear();
    int stat = 0;
    const off_t size = lseek(fd, 0, SEEK_END);
    if (size > 0) {
        off_t jmp = size - max_line_length;
        if (jmp < 0) jmp = 0;
        lseek(fd, jmp, SEEK_SET);

        char * raw = new char[max_line_length * sizeof(char)];
        const ssize_t n_read = read(fd, raw, max_line_length);
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
            buf.assign(chunk, begin, line_size - 1);
            if (! keep_flag) {
                stat = ftruncate(fd, size - line_size);
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
#else
    return 0;
#endif
}

int FileStack::pop_non_locked(string & buf) {
    DBG("");
    size_t size_after_pop = numeric_limits<size_t>::max();
    return pop_non_locked(buf, false, size_after_pop);
}

int FileStack::pop(string & buf, const bool keep_flag, size_t & size_after_pop) {
#ifndef WIN32
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
#else
    return 0;
#endif
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
        return -1;
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
        return -1;
    }
}

int FileStack::remove(const string & line) {
    DBG("");
#ifndef WIN32
    bool locked_internally = false;
    if (! locked) {
        lock();
        locked_internally = true;
    }
    const off_t size = lseek(fd, 0, SEEK_END);
    lseek(fd, 0, SEEK_SET);

    char * raw = new char[size * sizeof(char)];
    const ssize_t n_read = read(fd, raw, size);
    string buf;
    buf.assign(raw, n_read);
    delete [] raw;

    vector<string> tokens = split(buf, '\n');
    buf.clear();

    tokens.erase(std::remove(tokens.begin(), tokens.end(), line), tokens.end());

    lseek(fd, 0, SEEK_SET);
    int stat = ftruncate(fd, 0);
    for (auto it = tokens.begin(); it != tokens.end(); ++it) {
        buf = *it + '\n';
        size_t n = write(fd, buf.c_str(), buf.size());
    }

    if (locked_internally) {
        unlock();
    }
#endif
    return 0;
}

int64_t FileStack::push_non_locked(const string & buf) {
    DBG("");
    static const string nl("\n");
#ifndef WIN32
    lseek(fd, 0, SEEK_END);
    size_t n = write(fd, buf.c_str(), buf.size());
    if (buf.back() != nl.back()) {
        n += write(fd, nl.c_str(), nl.size());
    }
    return n;
#else
	return 0;
#endif
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

#ifndef WIN32
    size_t n_bytes, i;
    lseek(fd, 0, SEEK_SET);
    while ((n_bytes = read(fd, raw, chunk_size)) > 0) {
        for (i=0; i<n_bytes; i++) {
            if (raw[i] == '\n') {
                c++;
            }
        }
    }
#endif

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

#ifndef WIN32
    lseek(fd, 0, SEEK_SET);
    int stat = ftruncate(fd, 0);
#else
    int stat = 0;
#endif

    if (locked_internally) {
        unlock();
    }

    return stat;
}



bool FileStack::poll_query(const string & query, const double sleep_s, const size_t max_iter) {
    DBG("");
    string buf;
    const chrono::duration<double> sleep_time(sleep_s);
    for (size_t i=0; i < max_iter; ++i) {
        top(buf);
        if (buf.find(query) != string::npos) {
            DBG(string("") + " success : poll_iteration=" + to_string(i) + ", query=" + "\"" + query + "\"");
            return true;
        } else {
            DBG(string("") + " ongoing : poll_iteration=" + to_string(i) + ", query=" + "\"" + query + "\"");
        }
        if (buf.find("STOP") != string::npos) {
            throw(runtime_error("STOP on FileStack " + file_name));
        }
        this_thread::sleep_for(sleep_time);
    }
    throw(runtime_error("Could not discover keyword " + query + " on FileStack " + file_name
                      + " within " + to_string(double(max_iter) * sleep_s) + " seconds."));
    return false;  // TODO : finally decide on the semantics
};

bool FileStack::poll_size(const size_t size, const double sleep_s, const size_t max_iter) {
    DBG("");
    string buf;
    const chrono::duration<double> sleep_time(sleep_s);
    for (size_t i=0; i < max_iter; ++i) {
        if (this->size() == size) {
            return true;
        }
        this_thread::sleep_for(sleep_time);
    }
    throw(runtime_error("Could not detect size " + to_string(size) + " of FileStack " + file_name
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
