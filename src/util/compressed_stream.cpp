#include <stdexcept>
#include "compressed_stream.h"

Compressed_istream::Compressed_istream(const string & file_name) :
	file_name_(file_name),
	s_(file_name, std::ios_base::in | std::ios_base::binary)
{ }

size_t Compressed_istream::read(char * ptr, size_t count)
{
	s_.read(ptr, count);
	const size_t n = s_.gcount();
	if (n != count) {
		if (s_.eof())
			return n;
		else
			throw std::runtime_error("Error reading file " + file_name_);
	}
	return n;
}

void Compressed_istream::putback(char c)
{
	s_.putback(c);
	if (!s_.good())
		throw std::runtime_error("Error reading file " + file_name_);
}

Compressed_ostream::Compressed_ostream(const string & file_name, bool compressed) :
	file_name_(file_name),
	s_(file_name == "" ? auto_ptr<std::ostream>(&std::cout)
		: (compressed ? auto_ptr<std::ostream>(new zstr::ofstream(file_name)) : auto_ptr<std::ostream>(new std::ofstream(file_name.c_str()))))
{}

size_t Compressed_ostream::write(const char * ptr, size_t count)
{
	s_->write(ptr, count);
	if (!s_->good())
		throw std::runtime_error("Error writing file " + file_name_);
	return count;
}
