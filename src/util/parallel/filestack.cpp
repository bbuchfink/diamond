#include <string>
#include <iostream>
#include <chrono>
#include <thread>
#include <stdexcept>

#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "filestack.h"

using namespace std;


const string default_file_name = "default_stack.idx";
const int default_max_line_length = 4096;


FileStack::FileStack() : FileStack::FileStack(default_file_name) {
    DBG(__PRETTY_FUNCTION__);
    std::cerr << "FileStack: Using default file name " << default_file_name << std::endl;
}

FileStack::FileStack(const string & file_name) : FileStack::FileStack(file_name, default_max_line_length) {
    DBG(__PRETTY_FUNCTION__);
}

FileStack::FileStack(const string & file_name, int maximum_line_length) {
    DBG(__PRETTY_FUNCTION__);
    fd = open(file_name.c_str(), O_RDWR | O_CREAT, 00664);
    if (fd == -1) {
        throw(std::runtime_error("could not open file " + file_name));
    }
    this->locked = false;
    this->file_name = file_name;
    set_max_line_length(maximum_line_length);
}

FileStack::~FileStack() {
    DBG(__PRETTY_FUNCTION__);
    close(fd);
}

int FileStack::lock() {
    DBG(__PRETTY_FUNCTION__);
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
}

int FileStack::unlock() {
    DBG(__PRETTY_FUNCTION__);
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
}

int FileStack::pop(string & buf) {
    DBG(__PRETTY_FUNCTION__);
    return pop(buf, false);
}

int FileStack::top(string & buf) {
    DBG(__PRETTY_FUNCTION__);
    return pop(buf, true);
}

int FileStack::pop(int & i) {
    DBG(__PRETTY_FUNCTION__);
    string buf;
    const int get_status = pop(buf, false);
    if (get_status > 0) {
        return i = stoi(buf);
    } else {
        return -1;
    }
}

int FileStack::top(int & i) {
    DBG(__PRETTY_FUNCTION__);
    string buf;
    const int get_status = pop(buf, true);
    if (get_status > 0) {
        return i = stoi(buf);
    } else {
        return -1;
    }
}

int FileStack::pop(string & buf, const bool keep_flag) {
    DBG(__PRETTY_FUNCTION__);
    buf.clear();
    bool locked_internally = false;
    if (! locked) {
        lock();
        locked_internally = true;
    }
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
        // fprintf(stderr, "%d %d\n", begin, end);
        // fprintf(stderr, "chunk begin : %c|", chunk[begin]);
        // fprintf(stderr, "chunk end   : %c|", chunk[end]);
        // fprintf(stderr, "%s|", buf.c_str());
    }
    if (locked_internally) {
        unlock();
    }
    if (stat == -1) {
        return stat;
    } else {
        return buf.size();
    }
}

int FileStack::push(string buf) {
    DBG(__PRETTY_FUNCTION__);
    bool added_newline = false;
    if (buf.back() != '\n') {
        buf.append(1, '\n');
        added_newline = true;
    }
    bool locked_internally = false;
    if (! locked) {
        lock();
        locked_internally = true;
    }
    lseek(fd, 0, SEEK_END);
    const size_t n = write(fd, buf.c_str(), buf.size());
    if (locked_internally) {
        unlock();
    }
    if (added_newline) {
        buf.pop_back();
    }
    return n;
}

int FileStack::push(int i) {
    DBG(__PRETTY_FUNCTION__);
    string buf = to_string(i);
    return push(buf);
}

int FileStack::set_max_line_length(int n) {
    DBG(__PRETTY_FUNCTION__);
    const int minimum_line_length = 8;
    if (n < minimum_line_length) {
        n = minimum_line_length;
    }
    max_line_length = n;
    return max_line_length;
}

int FileStack::get_max_line_length() {
    DBG(__PRETTY_FUNCTION__);
    return max_line_length;
}

size_t FileStack::size() {
    DBG(__PRETTY_FUNCTION__);
    size_t n_bytes, i, c = 0;
    const size_t chunk_size = default_max_line_length;
    char * raw = new char[chunk_size * sizeof(char)];
    lock();
    lseek(fd, 0, SEEK_SET);
    while ((n_bytes = read(fd, raw, chunk_size)) > 0) {
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
    DBG(__PRETTY_FUNCTION__);
    lock();
    lseek(fd, 0, SEEK_SET);
    int stat = ftruncate(fd, 0);
    unlock();
    return stat;
}

bool FileStack::poll_query(const string & query, const double sleep_s, const size_t max_iter) {
    DBG(__PRETTY_FUNCTION__);
    string buf;
    const chrono::duration<double> sleep_time(sleep_s);
    for (size_t i=0; i < max_iter; ++i) {
        top(buf);
        if (buf.find(query) != string::npos) {
            DBG(string(__PRETTY_FUNCTION__) + " success : poll_iteration=" + to_string(i) + ", query=" + "\"" + query + "\"");
            return true;
        } else {
            DBG(string(__PRETTY_FUNCTION__) + " ongoing : poll_iteration=" + to_string(i) + ", query=" + "\"" + query + "\"");
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
    DBG(__PRETTY_FUNCTION__);
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
