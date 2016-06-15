#ifndef COMPRESSED_STREAM_H_
#define COMPRESSED_STREAM_H_

#include <string>
#include <memory>
#include "zstr.hpp"

using std::string;
using std::auto_ptr;

struct Compressed_istream
{
	Compressed_istream(const string &file_name);
	size_t read(char *ptr, size_t count);
	void putback(char c);
private:
	const string file_name_;
	zstr::ifstream s_;
};

struct Compressed_ostream
{
	Compressed_ostream(const string &file_name, bool compressed);
	size_t write(const char *ptr, size_t count);
	std::ostream& stream()
	{
		return *s_;
	}
	~Compressed_ostream();
private:
	const string file_name_;
	std::ostream *s_;
};

#endif