/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

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

#pragma once
#include <string>
#include <tuple>

#ifdef WIN32
const char PATH_SEPARATOR = '\\';
#else
const char PATH_SEPARATOR = '/';
#endif

#ifdef _WIN32
#define DEFAULT_LINE_DELIMITER "\r\n"
#else
#define DEFAULT_LINE_DELIMITER "\n"
#endif

enum class Color { RED, GREEN, YELLOW };

void set_color(Color color, bool err = false);
void reset_color(bool err = false);
std::string executable_path();
bool exists(const std::string &file_name);
void auto_append_extension(std::string &str, const char *ext);
std::string auto_append_extension_if_exists(const std::string &str, const char *ext);
size_t getCurrentRSS();
size_t getPeakRSS();
void log_rss();
size_t file_size(const char* name);
double total_ram();
std::tuple<char*, size_t, int> mmap_file(const char* filename);
void unmap_file(char* ptr, size_t size, int fd);
size_t l3_cache_size();
void mkdir(const std::string& dir);
void rmdir(const std::string& dir);
std::string containing_directory_absolute(const std::string& file_path);

#ifdef _MSC_VER
#define POPEN _popen
#define PCLOSE _pclose
#else
#define POPEN popen
#define PCLOSE pclose
#endif