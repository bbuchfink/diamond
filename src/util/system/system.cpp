#include <stdexcept>
#include <string.h>
#include "system.h"
#include "../string/string.h"

using namespace std;

#ifdef _MSC_VER
#include <windows.h>
#else
#include <unistd.h>
#include <sys/stat.h>
#endif

string executable_path() {
	char buf[4096];
#ifdef _MSC_VER
	if (GetModuleFileNameA(NULL, buf, sizeof(buf)) == 0)
		throw runtime_error("Error executing GetModuleFileNameA.");
	return string(buf);
#else
	if (readlink("/proc/self/exe", buf, sizeof(buf)) < 0)
		throw runtime_error("Error executing readlink on /proc/self/exe.");
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