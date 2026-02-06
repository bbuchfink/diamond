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

#pragma once

#include <map>
#include <mutex>
#include <string>
#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <unistd.h>
#include <fcntl.h>
#endif

#ifdef DEBUG
#define DBG(x) (std::cerr << __PRETTY_FUNCTION__ << ":" << __LINE__ << " " << x << std::endl)
#else
#define DBG(x)
#endif

class FileStack {
    public:
        FileStack();
        FileStack(const std::string & file_name);
        FileStack(const std::string & file_name, int maximum_line_length);
        ~FileStack();

        size_t size();

        int64_t pop(int64_t& i);
        int pop(std::string & buf);
        int pop(std::string & buf, size_t & size_after_pop);        

        int64_t top(int64_t& i);
        int top(std::string & buf);

        void remove(const std::string & line);

        int64_t push(int64_t i);
        int64_t push(const std::string & buf);
        int64_t push(const std::string & buf, size_t & size_after_push);
        int64_t fetch_add(int64_t n = 1);

        int get_max_line_length();
        int set_max_line_length(int n);

        int clear();
        
        int64_t seek(int64_t offset, int mode);
        size_t read(char* buf, size_t size);
        int64_t write(const char* buf, size_t size);
        int truncate(size_t size);

        bool poll_query(const std::string & query, const double sleep_s=0.5, const size_t max_iter=7200);
        bool poll_size(const size_t size, const double sleep_s=0.5, const size_t max_iter=7200);
        std::string file_name() const;
                
        int64_t pop_exclusive(std::string& buf);
        int64_t push_exclusive(const std::string& buf);

    private:
        
        int lock();
        int unlock();        
        
#ifdef _WIN32
        HANDLE hFile;
#else
        int fd;
        struct flock lck;
#endif
        std::string file_name_;
        off_t max_line_length;

        int pop(std::string &, const bool, size_t &);
        int64_t pop_non_locked(std::string &, const bool, size_t &);
        int64_t pop_non_locked(int64_t& i);
        int64_t push_non_locked(int64_t i);

        std::mutex mtx_;

};