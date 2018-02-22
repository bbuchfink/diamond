/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef SERIALIZER_H_
#define SERIALIZER_H_

#include <string>
#include "output_file.h"
#include "input_file.h"

using std::string;

struct Serializer
{

	Serializer(OutputFile &f):
		f_(f)
	{}

	Serializer& operator<<(int x)
	{
		f_.write(&x, 1);
		return *this;
	}

	Serializer& operator<<(const string &s)
	{
		f_.write(s.c_str(), s.length() + 1);
		return *this;
	}

	Serializer& operator<<(const vector<string> &v)
	{
		(*this) << (int)v.size();
		for (vector<string>::const_iterator i = v.begin(); i < v.end(); ++i)
			(*this) << *i;
		return *this;
	}

protected:

	OutputFile &f_;

};

struct Deserializer
{

	Deserializer():
		f_(NULL)
	{
	}

	Deserializer(InputFile &f) :
		f_(&f)
	{}

	Deserializer& operator>>(int &x)
	{
		if (f_->read(&x, 1) != 1)
			throw std::runtime_error("Unexpected end of buffer.");
		return *this;
	}

	Deserializer& operator>>(string &s)
	{
		if (!f_->read_until(s, '\0'))
			throw std::runtime_error("Unexpected end of buffer.");
		return *this;
	}

	Deserializer& operator>>(vector<string> &v)
	{
		int n;
		*this >> n;
		v.clear();
		string s;
		for (int i = 0; i < n; ++i) {
			*this >> s;
			v.push_back(s);
		}
		return *this;
	}

protected:

	InputFile *f_;

};

#endif