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

#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>
#include <errno.h>
#include "system.h"
#include "../string/string.h"
#include "../log_stream.h"
#ifdef _MSC_VER
  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
#else
  #include <unistd.h>
  #include <limits.h>
  #include <sys/stat.h>
  #include <sys/mman.h>
  #include <fcntl.h>
  #include <unistd.h>
  #ifdef __FreeBSD__
    #include <sys/types.h>
    #include <sys/sysctl.h>
  #elif defined(HAVE_SYSINFO)
    #include <sys/sysinfo.h>
  #endif
#endif

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::runtime_error;
using std::tuple;
using std::pair;
using std::make_pair;

string executable_path() {
	char buf[4096];
#ifdef _MSC_VER
	if (GetModuleFileNameA(NULL, buf, sizeof(buf)) == 0)
		throw runtime_error("Error executing GetModuleFileNameA.");
	return string(buf);
#else
	if (readlink("/proc/self/exe", buf, sizeof(buf)) < 0)
		throw std::runtime_error("Error executing readlink on /proc/self/exe.");
	return string(buf);
#endif
}

bool exists(const std::string &file_name) {
#ifdef _MSC_VER
	return GetFileAttributes(file_name.c_str()) != INVALID_FILE_ATTRIBUTES;
#else
	struct stat buffer;
	return stat(file_name.c_str(), &buffer) == 0;
#endif
}

size_t file_size(const char* name)
{
#ifdef WIN32
	HANDLE hFile = CreateFile(name, GENERIC_READ,
		FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, OPEN_EXISTING,
		FILE_ATTRIBUTE_NORMAL, NULL);
	if (hFile == INVALID_HANDLE_VALUE)
		return -1; // error condition, could call GetLastError to find out more

	LARGE_INTEGER size;
	if (!GetFileSizeEx(hFile, &size))
	{
		CloseHandle(hFile);
		return -1; // error condition, could call GetLastError to find out more
	}

	CloseHandle(hFile);
	return size.QuadPart;
#else
	struct stat stat_buf;
	int rc = stat(name, &stat_buf);
	return rc == 0 ? stat_buf.st_size : -1;
#endif
}

void auto_append_extension(string &str, const char *ext)
{
	if (!str.empty() && !ends_with(str, ext))
		str += ext;
}

string auto_append_extension_if_exists(const string &str, const char *ext) {
	if (!ends_with(str, ext) && exists(str + ext))
		return str + ext;
	return str;
}

void log_rss() {
	log_stream << "Current RSS: " << convert_size(getCurrentRSS()) << ", Peak RSS: " << convert_size(getPeakRSS()) << endl;
}

void set_color(Color color, bool err) {
#ifdef WIN32
	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	WORD c = FOREGROUND_RED | FOREGROUND_BLUE | FOREGROUND_GREEN;
	switch (color) {
	case Color::RED:
		c = FOREGROUND_RED;
		break;
	case Color::GREEN:
		c = FOREGROUND_GREEN;
		break;
	case Color::YELLOW:
		c = FOREGROUND_RED | FOREGROUND_GREEN;
		break;
	default:
		break;
	}
	SetConsoleTextAttribute(hConsole, c);
#else
	auto &s = err ? cerr : cout;
	s << "\033[";
	switch (color) {
	case Color::RED:
		s << 31;
		break;
	case Color::GREEN:
		s << 32;
		break;
	case Color::YELLOW:
		s << "1;33";
		break;
	default:
		break;
	}
	s << "m";
#endif
}

void reset_color(bool err) {
#ifdef WIN32
	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_BLUE | FOREGROUND_GREEN);
#else
	(err ? cerr : cout) << "\033[" << "0;39" << 'm';
#endif
}

double total_ram() {
#ifdef __FreeBSD__
	int mib[2] = { CTL_HW, HW_REALMEM };
	u_int namelen = sizeof(mib) / sizeof(mib[0]);
	uint64_t oldp;
	size_t oldlenp = sizeof(oldp);

	if (sysctl(mib, namelen, &oldp, &oldlenp, NULL, 0) < 0)
		return 0.0;
	else
		return oldp / 1e9;
#elif defined(HAVE_SYSINFO)
	struct sysinfo info;
	if (sysinfo(&info) != 0)
		return 0.0;
	return (double)info.totalram / 1e9;
#else
	return 0.0;
#endif
}

tuple<char*, size_t, int> mmap_file(const char* filename) {
#ifdef WIN32
	throw std::runtime_error("Memory mapping not supported on Windows.");
#else
	void* addr;
	int fd;
	struct stat sb;
	size_t length;

	fd = open(filename, O_RDONLY);
	if (fd == -1)
		throw std::runtime_error(string("Error opening file: ") + filename);

	if (fstat(fd, &sb) == -1)
		throw std::runtime_error(string("Error calling fstat on file: ") + filename);

	length = sb.st_size;

	addr = mmap(NULL, length, PROT_READ, MAP_SHARED, fd, 0);
	if (addr == MAP_FAILED)
		throw std::runtime_error(string("Error calling mmap on file: ") + filename);

	return std::tuple<char*, size_t, int>((char*)addr, length, fd);
#endif
}

void unmap_file(char* ptr, size_t size, int fd) {
#ifdef WIN32
#else
	munmap((void*)ptr, size);
	close(fd);
#endif
}

size_t l3_cache_size() {
#if defined(_MSC_VER) || defined(__APPLE__) || defined(__FreeBSD__) || (! defined(_SC_LEVEL3_CACHE_SIZE))
	return 0;
#else
	const auto s = sysconf(_SC_LEVEL3_CACHE_SIZE);
	return s == -1 ? 0 : s;
#endif
}

void mkdir(const std::string& dir) {
#ifdef WIN32
	CreateDirectory(dir.c_str(), NULL);
#else
	errno = 0;
	if (::mkdir(dir.c_str(), 493) != 0) {
		if (errno == EEXIST) {
		}
		else
			throw(std::runtime_error("could not create temporary directory " + dir));
	}
#endif
}

void rmdir(const std::string& path)
{
#ifdef _WIN32
	const int wlen = MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, nullptr, 0);
	if (wlen <= 0) throw runtime_error("MultiByteToWideChar");

	std::wstring wpath(wlen, L'\0');
	MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, &wpath[0], wlen);

	RemoveDirectoryW(wpath.c_str());
#else
	::rmdir(path.c_str());
#endif
}

static bool is_sep_char(char c) {
#ifdef _WIN32
	return c == '\\' || c == '/';
#else
	return c == '/';
#endif
}

static bool ends_with_sep(const std::string& s) {
	return !s.empty() && is_sep_char(s.back());
}

static std::string last_component(const std::string& p) {
	if (p.empty()) return std::string();
	std::size_t end = p.size();
	while (end > 0 && is_sep_char(p[end - 1])) --end;
	std::size_t start = end;
	while (start > 0 && !is_sep_char(p[start - 1])) --start;
	return p.substr(start, end - start);
}

#ifndef _WIN32

static bool is_abs_posix(const std::string& p) {
	return !p.empty() && p[0] == '/';
}

static std::string get_cwd_posix() {
	std::vector<char> buf(256);
	for (;;) {
		if (getcwd(buf.data(), buf.size())) return std::string(buf.data());
		buf.resize(buf.size() * 2);
	}
}

static std::string lex_normalize_posix(const std::string& path) {
	const bool abs = is_abs_posix(path);
	std::vector<std::string> stack;
	std::size_t i = 0, n = path.size();

	while (i < n) {
		while (i < n && path[i] == '/') ++i;
		std::size_t start = i;
		while (i < n && path[i] != '/') ++i;
		std::string token = path.substr(start, i - start);
		if (token.empty() || token == ".") continue;
		if (token == "..") {
			if (!stack.empty()) stack.pop_back();
			else if (!abs) stack.push_back("..");
		}
		else {
			stack.push_back(token);
		}
	}

	std::string out;
	if (abs) out.push_back('/');
	for (std::size_t k = 0; k < stack.size(); ++k) {
		if (k) out.push_back('/');
		out += stack[k];
	}
	if (out.empty()) return abs ? "/" : ".";
	return out;
}

static std::string parent_dir_posix(const std::string& abs_path) {
	if (abs_path == "/") return "/";
	std::size_t pos = abs_path.find_last_of('/');
	if (pos == std::string::npos) return ".";
	if (pos == 0) return "/";
	return abs_path.substr(0, pos);
}

#else

static std::wstring widen_utf8(const std::string& s) {
	if (s.empty()) return std::wstring();
	int wlen = MultiByteToWideChar(CP_UTF8, 0, s.data(), (int)s.size(), nullptr, 0);
	std::wstring w(wlen, L'\0');
	MultiByteToWideChar(CP_UTF8, 0, s.data(), (int)s.size(), &w[0], wlen);
	return w;
}

static std::string narrow_utf8(const std::wstring& w) {
	if (w.empty()) return std::string();
	int len = WideCharToMultiByte(CP_UTF8, 0, w.data(), (int)w.size(), nullptr, 0, nullptr, nullptr);
	std::string s(len, '\0');
	WideCharToMultiByte(CP_UTF8, 0, w.data(), (int)w.size(), &s[0], len, nullptr, nullptr);
	return s;
}

#endif

pair<string, string> absolute_path(const std::string& file_path) {
	const std::string fp = file_path.empty() ? "." : file_path;
	const std::string base = last_component(fp);
	const bool treat_as_dir = ends_with_sep(fp) || base == "." || base == "..";

#ifndef WIN32
	std::string joined;
	if (is_abs_posix(fp)) {
		joined = fp;
	}
	else {
		std::string cwd = get_cwd_posix();
		if (cwd.empty()) return {};
		joined = cwd + "/" + fp;
	}

	std::string abs_norm = lex_normalize_posix(joined);
	if (treat_as_dir) return make_pair(abs_norm, string());
	return make_pair(parent_dir_posix(abs_norm), base);
#else
	std::wstring winput = widen_utf8(fp);
	if (winput.empty()) winput = L".";

	DWORD need = GetFullPathNameW(winput.c_str(), 0, nullptr, nullptr);
	if (need == 0) return {};

	std::vector<wchar_t> buf(need);
	DWORD written = GetFullPathNameW(winput.c_str(), (DWORD)buf.size(), buf.data(), nullptr);
	if (written == 0) return {};

	std::wstring full(buf.data(), written);

	for (auto& ch : full) if (ch == L'/') ch = L'\\';

	if (treat_as_dir) {
		return make_pair(narrow_utf8(full), string());
	}

	std::size_t pos = full.find_last_of(L"\\/");
	if (pos == std::wstring::npos) return make_pair(narrow_utf8(full), base);

	if (pos == 2 && full.size() >= 3 && full[1] == L':') ++pos;
	std::wstring parent = full.substr(0, pos);
	return make_pair(narrow_utf8(parent), base);
#endif
}

bool stdout_is_a_tty() {
#if defined(_WIN32)
	HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
	if (h == nullptr || h == INVALID_HANDLE_VALUE) return false;
	if (GetFileType(h) != FILE_TYPE_CHAR) return false;
	DWORD mode;
	return GetConsoleMode(h, &mode) != 0;
#else
	return ::isatty(STDOUT_FILENO) == 1;
#endif
}

bool is_absolute_path(const std::string& path) {
	if (path.empty())
		return false;
	if (path[0] == '/')
		return true;
	if (path[0] == '\\')
		return true;
	if (path.size() >= 3 &&
		std::isalpha(static_cast<unsigned char>(path[0])) &&
		path[1] == ':' &&
		(path[2] == '/' || path[2] == '\\'))
	{
		return true;
	}
	return false;
}