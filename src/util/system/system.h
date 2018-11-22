#ifndef UTIL_SYSTEM_SYSTEM_H_
#define UTIL_SYSTEM_SYSTEM_H_

#include <string>

std::string executable_path();
bool exists(const std::string &file_name);
void auto_append_extension(std::string &str, const char *ext);
void auto_append_extension_if_exists(std::string &str, const char *ext);

#endif