/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
void mkdir(const std::string& dir);
void rmdir(const std::string& dir);
std::pair<std::string, std::string> absolute_path(const std::string& file_path);
bool is_absolute_path(const std::string& path);
void remove_tmp_file(const std::string& file_name);
std::string create_temp_directory(const std::string& parent_dir, const std::string& prefix, unsigned max_attempts = 100);
long max_open_files_per_process();
long raise_open_files_limit(long desired = 0);

inline std::string containing_directory(const std::string& file_name) {
	return file_name.substr(0, file_name.find_last_of(PATH_SEPARATOR));
}


#ifdef _MSC_VER
#define POPEN _popen
#define PCLOSE _pclose
#else
#define POPEN popen
#define PCLOSE pclose
#endif