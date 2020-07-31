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
#include <sys/sysinfo.h>
#endif

using std::string;
using std::cout;
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

void set_color(Color color) {
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
	cout << "\033[";
	switch (color) {
	case Color::RED:
		cout << 31;
		break;
	case Color::GREEN:
		cout << 32;
		break;
	default:
		break;
	}
	cout << 'm';
#endif
}

void reset_color() {
#ifdef WIN32
	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_BLUE | FOREGROUND_GREEN);
#else
	cout << "\033[" << 39 << 'm';
#endif
}

double total_ram() {
#ifdef WIN32
	return 0.0;
#else
	sysinfo info;
	if (sysinfo(&info) != 0)
		return 0.0;
	return (double)info.totalram / double(1 << 30);
#endif
}