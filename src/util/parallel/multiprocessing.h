/****
DIAMOND protein aligner
Copyright (C) 2019-2024 Max Planck Society for the Advancement of Science e.V.

Code developed by Klaus Reuter

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#pragma once
// various tools related to the multiprocessing parallelization, to be moved elsewhere

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>

std::vector<std::string> split(const std::string & str, const char delim);
std::string join(const std::vector<std::string> & tokens, const char delim);
std::string quote(const std::string & str);
std::string unquote(const std::string & str);
void copy(const std::string & src_file_name, const std::string & dst_file_name);
std::string join_path(const std::string & path_1, const std::string & path_2);
bool file_exists(const std::string & file_name);


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
