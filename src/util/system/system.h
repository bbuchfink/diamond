#pragma once
#include <stdio.h>
#include <string>
#include <tuple>

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

#ifdef _MSC_VER
#define POPEN _popen
#define PCLOSE _pclose
#else
#define POPEN popen
#define PCLOSE pclose
#endif