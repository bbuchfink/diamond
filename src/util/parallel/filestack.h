/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink, Klaus Reuter
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
// SPDX-License-Identifier: GPL-3.0-or-later

#pragma once

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
        void close();
        template<typename Job>
        FileStack(const std::string& file_name, Job& job):
            FileStack(file_name)
        {
            job.register_sync_file(file_name);
		}

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
        void remove();

    private:
        
        int lock();
        int unlock();
        size_t size_unlocked();
        
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