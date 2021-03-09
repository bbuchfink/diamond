#ifndef _MULTIPROCESSING_H_
#define _MULTIPROCESSING_H_
// various tools related to the multiprocessing parallelization, to be moved elsewhere

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>

using std::vector;
using std::string;

vector<string> split(const string & str, const char delim);
string join(const vector<string> & tokens, const char delim);
string quote(const string & str);
string unquote(const string & str);
void copy(const string & src_file_name, const string & dst_file_name);
string join_path(const string & path_1, const string & path_2);
bool file_exists(const string & file_name);


template<typename T>
std::istream& load_scalar(std::istream& is, T & v)
{
    is.read(reinterpret_cast<char*>(&v), sizeof(v));
    return is;
}

template<typename S>
std::istream& load_string(std::istream& is, S & s)
{
    auto size = s.size();
    is.read(reinterpret_cast<char*>(&size), sizeof(size));
	std::vector<char> v;
    v.resize(size);
    is.read(reinterpret_cast<char*>(v.data()), size * sizeof(char));
	v.push_back('\0');
	s = reinterpret_cast<char*>(v.data());
    return is;
}

template<typename CONTAINER>
std::istream& load_vector(std::istream& is, CONTAINER & v)
{
    auto size = v.size();
    is.read(reinterpret_cast<char*>(&size), sizeof(size));
    v.resize(size);
    is.read(reinterpret_cast<char*>(v.data()), v.size() * sizeof(typename CONTAINER::value_type));
    return is;
}

template<typename T>
std::ostream& save_scalar(std::ostream& os, T & v)
{
    os.write(reinterpret_cast<char const*>(&v), sizeof(v));
    return os;
}

template<typename CONTAINER>
std::ostream& save_vector(std::ostream& os, CONTAINER const& v)
{
    auto size = v.size();
    os.write(reinterpret_cast<char const*>(&size), sizeof(size));
    os.write(reinterpret_cast<char const*>(v.data()), v.size() * sizeof(typename CONTAINER::value_type));
    return os;
}

template<typename T>
std::string append_label(const std::string & str, const T & label, const size_t width=6) {
	std::stringstream ss;
	ss << std::setfill('0') << std::setw(width) << label;
	return str + ss.str();
}

#endif
