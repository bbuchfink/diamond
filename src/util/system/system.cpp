#include <stdexcept>
#include <string.h>
#include "system.h"

using namespace std;

#ifdef _MSC_VER
#include <windows.h>
#else
#include <unistd.h>
#include <sys/stat.h>
#include <linux/limits.h>
#endif

string executable_path() {
#ifdef _MSC_VER
	char buf[4096];
	if (GetModuleFileNameA(NULL, buf, sizeof(buf)) == 0)
		throw runtime_error("Error executing GetModuleFileNameA.");
	return string(buf);
#else
	char buf[PATH_MAX + 1];
	if (readlink("/proc/self/exe", buf, PATH_MAX) < 0)
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
	size_t l = strlen(ext);
	if (str.length() < l || (str.length() >= l && str.substr(str.length() - l, string::npos) != ext))
		str += ext;
}

void auto_append_extension_if_exists(string &str, const char *ext) {

}