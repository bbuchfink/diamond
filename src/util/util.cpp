#include <stdio.h>
#include <exception>
#include <sstream>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include "../basic/config.h"
#include "log_stream.h"
#include "complexity_filter.h"
#include "util.h"
#include "temp_file.h"
#include "binary_file.h"

Message_stream message_stream;
Message_stream verbose_stream (false);
Message_stream log_stream (false);

const Complexity_filter Complexity_filter::instance;

#ifndef _MSC_VER
const char dir_separator = '/';
#else
const char dir_separator = '\\';
#endif

string extract_dir(const string & s)
{
	return s.find_last_of(dir_separator) == string::npos ? "" : s.substr(0, s.find_last_of(dir_separator));
}

unsigned Temp_file::n = 0;
unsigned long Temp_file::hash_key;

Temp_file::Temp_file()
{
	if (n == 0) {
#ifdef WIN32
		LARGE_INTEGER count;
		QueryPerformanceCounter(&count);
		hash_key = (unsigned long)(count.HighPart + count.LowPart + count.QuadPart + GetCurrentProcessId());
#else
		timeval count;
		gettimeofday(&count, NULL);
		hash_key = count.tv_sec + count.tv_usec + getpid();
#endif
	}
	std::stringstream ss;
	ss.setf(std::ios::hex, std::ios::basefield);
	if (config.tmpdir != "")
		ss << config.tmpdir << dir_separator;
	ss << "diamond-" << hash_key << "-" << n++ << ".tmp";
	ss >> this->file_name_;
	this->f_ = fopen(this->file_name_.c_str(), "w+b");
	if (this->f_ == 0)
		throw std::runtime_error("Error opening temporary file: " + this->file_name_);
#ifndef _MSC_VER
	unlink(this->file_name_.c_str());
#endif
}

string Temp_file::get_temp_dir()
{
	Temp_file t;
	Input_stream f(t);
	f.close_and_delete();
	return extract_dir(f.file_name);
}