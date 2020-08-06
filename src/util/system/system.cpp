#include <stdexcept>
#include <string.h>
#include <iostream>
#include "system.h"
#include "../string/string.h"
#include "../log_stream.h"

#ifdef _MSC_VER
  #include <windows.h>
#else
  #include <unistd.h>
  #include <sys/stat.h>
  #ifndef  __APPLE__
    #ifdef __FreeBSD__
      #include <sys/types.h>
      #include <sys/sysctl.h>
    #else
      #include <sys/sysinfo.h>
    #endif
  #endif
#endif

using std::string;
using std::cout;
using std::cerr;
using std::endl;

string executable_path() {
	char buf[4096];
#ifdef _MSC_VER
	if (GetModuleFileNameA(NULL, buf, sizeof(buf)) == 0)
		throw std::runtime_error("Error executing GetModuleFileNameA.");
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
	if (!ends_with(str, ext))
		str += ext;
}

void auto_append_extension_if_exists(string &str, const char *ext) {
	if (!ends_with(str, ext))
		if (exists(str + ext))
			str += ext;
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
#if defined(WIN32) || defined(__APPLE__)
	return 0.0;
#elif defined(__FreeBSD__)
	int mib[2] = { CTL_HW, HW_REALMEM };
	u_int namelen = sizeof(mib) / sizeof(mib[0]);
	uint64_t oldp;
	size_t oldlenp = sizeof(oldp);

	if (sysctl(mib, namelen, &oldp, &oldlenp, NULL, 0) < 0)
		return 0.0;
	else
		return oldp / 1e9;
#else
	struct sysinfo info;
	if (sysinfo(&info) != 0)
		return 0.0;
	return (double)info.totalram / 1e9;
#endif
}