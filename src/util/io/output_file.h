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

#ifndef OUTPUT_FILE_H_
#define OUTPUT_FILE_H_

#include <map>
#include <stdexcept>
#include <string>
#include "sink.h"
#include "file_sink.h"
#include "compressed_stream.h"
#include "../text_buffer.h"

using std::string;

struct OutputFile
{
	OutputFile()
	{}
	OutputFile(string &file_name, bool compressed = false):
		file_name_(file_name),
		sink_(new FileSink(file_name))
	{
		if (compressed)
			sink_ = new ZlibSink(sink_);
	}
	void close()
	{
		sink_->close();
	}
	template<typename _t>
	void write(const _t *ptr, size_t count)
	{
		sink_->write((const char*)ptr, sizeof(_t)*count);
	}
	template<class _t>
	void write(const vector<_t> &v)
	{
		write(v.data(), v.size());
	}
	void write_c_str(const string &s);
	void seekp(size_t p);
	void remove();
	size_t tell();

	template<typename _k, typename _v>
	void write_map_csv(typename std::map<_k, _v>::const_iterator begin, typename std::map<_k, _v>::const_iterator end)
	{
		TextBuffer buf;
		for (typename std::map<_k, _v>::const_iterator i = begin; i != end; ++i) {
			buf << i->first << '\t' << i->second << '\n';
			write(buf.get_begin(), buf.size());
			buf.clear();
		}

	}

	~OutputFile()
	{
		delete sink_;
	}
	const string file_name() const
	{
		return file_name_;
	}
	FILE* file() const
	{
		return sink_->file();
	}
protected:
	string file_name_;
	Sink *sink_;
};

#endif