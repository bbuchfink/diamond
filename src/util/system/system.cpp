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

#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>
#include <errno.h>
#include <cstring>
#include <system_error>
#include <random>
#include <cstdint>

#include "system.h"
#include "../string/string.h"
#include "../log_stream.h"
#include "basic/config.h"

#ifdef _WIN32
  #define NOMINMAX
  #define WIN32_LEAN_AND_MEAN
  #include <windows.h>
  #include <cstdio>          // _getmaxstdio
#else
  #include <unistd.h>
  #include <limits.h>
  #include <sys/stat.h>
  #include <sys/mman.h>
  #include <fcntl.h>
  #include <unistd.h>
  #include <sys/resource.h>  // getrlimit, RLIMIT_NOFILE
  #if defined(__APPLE__)
     #include <sys/sysctl.h>  // sysctlbyname (kern.maxfilesperproc)
  #endif
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

long max_open_files_per_process() {
#if defined(_WIN32)	
	return static_cast<long>(_getmaxstdio());
#else
	struct rlimit rl;
	if (getrlimit(RLIMIT_NOFILE, &rl) != 0)
		return -1;
	if (rl.rlim_cur == RLIM_INFINITY)
		return -1;
	return static_cast<long>(rl.rlim_cur);
#endif
}

long raise_open_files_limit(long desired) {
#if defined(_WIN32)
	int target = (desired > 0) ? static_cast<int>(desired) : 8192;
	if (target > 8192) target = 8192;
	if (_setmaxstdio(target) == -1)
		return -1;
	return static_cast<long>(_getmaxstdio());
#else
	struct rlimit rl;
	if (getrlimit(RLIMIT_NOFILE, &rl) != 0)
		return -1;

	rlim_t target = (desired > 0) ? static_cast<rlim_t>(desired) : rl.rlim_max;	
	if (rl.rlim_max != RLIM_INFINITY && target > rl.rlim_max)
		target = rl.rlim_max;

#if defined(__APPLE__)	
	int maxproc = 0;
	size_t sz = sizeof(maxproc);
	if (sysctlbyname("kern.maxfilesperproc", &maxproc, &sz, nullptr, 0) == 0
		&& maxproc > 0 && target > static_cast<rlim_t>(maxproc))
		target = static_cast<rlim_t>(maxproc);
#endif

	rl.rlim_cur = target;
	if (setrlimit(RLIMIT_NOFILE, &rl) != 0)
		return -1;

	if (getrlimit(RLIMIT_NOFILE, &rl) != 0)
		return -1;
	return static_cast<long>(rl.rlim_cur);
#endif
}

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
		throw runtime_error(string("Error getting size on file: ") + name + ": " + std::system_category().message(GetLastError()));

	LARGE_INTEGER size;
	if (!GetFileSizeEx(hFile, &size))
	{
		CloseHandle(hFile);
		throw runtime_error(string("Error getting size on file: ") + name + ": " + std::system_category().message(GetLastError()));
	}

	CloseHandle(hFile);
	return size.QuadPart;
#else
	struct stat stat_buf;
	int rc = stat(name, &stat_buf);
	if (rc != 0)
		throw runtime_error(string("Error getting size on file: ") + name + ": " + std::strerror(errno));
	return stat_buf.st_size;
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
	*log_stream << "Current RSS: " << convert_size(getCurrentRSS()) << ", Peak RSS: " << convert_size(getPeakRSS()) << endl;
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

void remove_tmp_file(const std::string& file_name) {
	if (config.keep_temp_files)
		return;
	if(remove(file_name.c_str()) != 0)
		fprintf(stderr, "Warning: Failed to delete temporary file %s\n", file_name.c_str());
}

static std::string join_path(const std::string& dir, const std::string& name) {
	if (dir.empty()) {
		return name;
	}

	char last = dir[dir.size() - 1];

#ifdef _WIN32
	if (last == '/' || last == '\\') {
		return dir + name;
	}
	return dir + "\\" + name;
#else
	if (last == '/') {
		return dir + name;
	}
	return dir + "/" + name;
#endif
}

static std::string random_suffix(std::mt19937_64& rng) {
	static const char alphabet[] =
		"0123456789"
		"ABCDEFGHIJKLMNOPQRSTUVWXYZ"
		"abcdefghijklmnopqrstuvwxyz";

	std::uniform_int_distribution<std::size_t> dist(
		0, sizeof(alphabet) - 2
	);

	std::string s;
	s.reserve(12);

	for (int i = 0; i < 12; ++i) {
		s.push_back(alphabet[dist(rng)]);
	}

	return s;
}

#ifdef _WIN32

static std::wstring utf8_to_wide(const std::string& s) {
	if (s.empty()) {
		return std::wstring();
	}
	int needed = MultiByteToWideChar(
		CP_UTF8,
		MB_ERR_INVALID_CHARS,
		s.data(),
		static_cast<int>(s.size()),
		NULL,
		0
	);
	if (needed <= 0) {
		throw std::runtime_error("invalid UTF-8 path");
	}
	std::wstring out(static_cast<std::size_t>(needed), L'\0');
	int written = MultiByteToWideChar(
		CP_UTF8,
		MB_ERR_INVALID_CHARS,
		s.data(),
		static_cast<int>(s.size()),
		&out[0],
		needed
	);
	if (written != needed) {
		throw std::runtime_error("UTF-8 to UTF-16 conversion failed");
	}
	return out;
}

#endif

std::string create_temp_directory(const std::string& parent_dir, const std::string& prefix, unsigned max_attempts) {
	if (parent_dir.find('\0') != std::string::npos) {
		throw std::invalid_argument("parent_dir contains NUL byte");
	}	
	std::random_device rd;
	std::mt19937_64 rng((static_cast<std::uint64_t>(rd()) << 32) ^ static_cast<std::uint64_t>(rd()));
	for (unsigned attempt = 0; attempt < max_attempts; ++attempt) {
		std::string name = prefix + random_suffix(rng);
		std::string path = join_path(parent_dir, name);
#ifdef _WIN32
		std::wstring wide_path = utf8_to_wide(path);
		if (CreateDirectoryW(wide_path.c_str(), NULL)) {
			return path;
		}
		DWORD err = GetLastError();
		if (err == ERROR_ALREADY_EXISTS || err == ERROR_FILE_EXISTS) {
			continue;
		}
		throw std::system_error(static_cast<int>(err), std::system_category(), "CreateDirectoryW failed for " + path);
#else
		if (::mkdir(path.c_str(), 0700) == 0) {
			return path;
		}
		int err = errno;
		if (err == EEXIST) {
			continue;
		}
		throw std::system_error(err, std::generic_category(), "mkdir failed for " + path);
#endif
	}
	throw runtime_error("could not create unique temporary directory after many attempts");
}