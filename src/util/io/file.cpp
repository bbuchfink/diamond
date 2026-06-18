/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#include <errno.h>
#ifdef _WIN32
#include <io.h>
#include <fcntl.h>
#endif
#include <string.h> // strerror
#include "file.h"
#include "temp_file.h"
#include "util/data_structures/mem_buffer.h"
#include "util/system/system.h"
#include "util/enum.h"
#include "basic/config.h"

using std::string;
using std::runtime_error;

static CompressionLib detect_compressor(const char* b) {
	if ((b[0] == '\x1F' && b[1] == '\x8B')         // gzip header
		|| (b[0] == '\x78' && (b[1] == '\x01'      // zlib header
			|| b[1] == '\x9C'
			|| b[1] == '\xDA')))
		return CompressionLib::ZLIB;
	if (b[0] == '\x28' && b[1] == '\xb5' && b[2] == '\x2f' && b[3] == '\xfd')
		return CompressionLib::ZSTD;
	return CompressionLib::NONE;
}

File::File(Temporary):
	line_buf_((char*)malloc(1024)),
	line_buf_size_(1024),
	decompressor_(new PassThrough),
	compressor_(new PassThroughCompressor)
{	
	TempFileData d = TempFile::init(true);
#ifdef _MSC_VER
	file_ = fopen(d.name.c_str(), "w+b");
#else
	file_ = fdopen(d.fd, "w+b");
#endif
	if (file_ == 0)
		throw runtime_error("Error opening temporary file " + d.name + ". " + strerror(errno));
	if (setvbuf(file_, nullptr, _IOFBF, config.file_buffer_size) != 0)
		throw runtime_error("Error setting buffer size for temporary file " + d.name + ". " + strerror(errno));
	unlinked_ = d.unlinked;
	file_name_ = d.name;
	auto_delete_ = true;
	seekable_ = true;
}

File::File(const string& name, const char* mode, Flags flags, CompressionLib compression):
	line_buf_((char*)malloc(1024)),
	line_buf_size_(1024),
	decompressor_(new PassThrough),
	compressor_(new PassThroughCompressor)
{
	if (name.empty()) {
		if (flag_any(flags, Flags::TREAT_BLANK_AS_STDIN)) {
			file_ = stdin;
#ifdef _WIN32
			_setmode(_fileno(stdin), _O_BINARY);
#endif
		}
		else if (flag_any(flags, Flags::TREAT_BLANK_AS_STDOUT)) {
			file_ = stdout;
#ifdef _WIN32
			_setmode(_fileno(stdout), _O_BINARY);
#endif
		}
	}
	else
		file_ = fopen(name.c_str(), mode);
	if (file_ == 0)
		throw runtime_error("Error opening file " + name + ". " + strerror(errno));
	if (setvbuf(file_, nullptr, _IOFBF, config.file_buffer_size) != 0)
		throw runtime_error("Error setting buffer size for file " + name + ". " + strerror(errno));
	unlinked_ = false;
	file_name_ = name;
	auto_delete_ = false;
	seekable_ = true;
	try {
		seek(1, SEEK_CUR);
		seek(-1, SEEK_CUR);
	}
	catch (const runtime_error&) {
		seekable_ = false;
	}
	if (strcmp(mode, "rb") == 0) {
		if (flag_any(flags, Flags::DETECT_COMPRESSION)) {
			const string magic = peek(4);
			if (magic.length() == 4) {
				switch (detect_compressor(magic.c_str())) {
				case CompressionLib::ZLIB:
					decompressor_.reset(new ZlibDecompressor);
					break;
				case CompressionLib::ZSTD:
#ifdef WITH_ZSTD
					decompressor_.reset(new ZstdDecompressor);
#endif
					break;
				default:
					break;
				}
			}
		}
	}
	else if (strcmp(mode, "wb") == 0 && compression != CompressionLib::NONE) {
		switch (compression) {
		case CompressionLib::ZLIB:
			compressor_.reset(new ZlibCompressor);
			break;
		case CompressionLib::ZSTD:
#ifdef WITH_ZSTD
			compressor_.reset(new ZstdCompressor);
#endif
			break;
		default:
			;
		}
	}
}

File::File(File&& f) noexcept {
	file_ = f.file_;
	unlinked_ = f.unlinked_;
	file_name_ = f.file_name_;
	auto_delete_ = f.auto_delete_;
	line_buf_ = f.line_buf_;
	line_buf_size_ = f.line_buf_size_;
	decompressor_ = std::move(f.decompressor_);
	compressor_ = std::move(f.compressor_);
	f.file_ = nullptr;
	f.line_buf_ = nullptr;
	f.line_buf_size_ = 0;
}

File& File::operator=(File&& f) noexcept {
	close();
	free(line_buf_);
	file_ = f.file_;
	unlinked_ = f.unlinked_;
	file_name_ = f.file_name_;
	auto_delete_ = f.auto_delete_;
	line_buf_ = f.line_buf_;
	line_buf_size_ = f.line_buf_size_;
	decompressor_ = std::move(f.decompressor_);
	compressor_ = std::move(f.compressor_);
	f.file_ = nullptr;
	f.line_buf_ = nullptr;
	f.line_buf_size_ = 0;
	return *this;
}

void File::close() {
	if (compressor_)
		compressor_->close(file_);
	if (file_) {
		fclose(file_);
		if (auto_delete_ && !unlinked_)
			remove_tmp_file(file_name_);
	}
	file_ = nullptr;
}

File::~File() {
	free(line_buf_);
	close();
}

void File::write(const void* ptr, size_t n) {
	if (compressor_->fwrite(ptr, 1, n, file_) != n) {
		perror(0);
		throw std::runtime_error("Error writing to file " + file_name_);
	}
}

void File::write_c_str(const char* s) {
	write(s, strlen(s) + 1);
}

void File::seek(int64_t p, int origin)
{
	if(decompressor_ && decompressor_->lib() != CompressionLib::NONE && (p != 0 || origin != SEEK_SET))
		throw runtime_error("Seeking is not supported for compressed files.");
#ifdef WIN32
	if (_fseeki64(file_, p, origin) != 0) {
		perror(0);
		throw std::runtime_error("Error calling fseek.");
	}
#else
	if (fseek(file_, p, origin) != 0) {
		perror(0);
		throw std::runtime_error("Error calling fseek.");
	}
#endif
}

void File::seek(int64_t p)
{
	seek(p, SEEK_SET);
}

void File::rewind()
{
	seek(0, SEEK_SET);
	clearerr(file_);
	line_buf_[0] = 0;
	line_count = 0;
	decompressor_->reset();
}

int64_t File::tell()
{	
#ifdef WIN32
	int64_t x;
	if ((x = _ftelli64(file_)) == (int64_t)-1)
		throw std::runtime_error("Error executing ftell on stream " + file_name_);
	return x;
#else
	const long n = ftell(file_);
	if (n < 0) {
		perror(0);
		throw std::runtime_error("Error calling ftell.");
	}
	return n;
#endif
}

int64_t File::size() {
	const int64_t pos = tell();
	seek(0l, SEEK_END);
	const int64_t s = tell();
	seek(pos, SEEK_SET);
	return s;
}

FILE* File::file() {
	return file_;
}

void File::read(void* ptr, size_t n) {
	const size_t r = decompressor_->fread(ptr, 1, n, file_);
	if(r != n)
		throw runtime_error("Error reading file " + file_name_ + ". " + strerror(errno));
}

size_t File::read_max(void* ptr, size_t n) {
	return decompressor_->fread(ptr, 1, n, file_);
}

const char* File::read_bytes(size_t n) {
	static MemBuffer<char> buf;
	buf.resize(n);
	read(buf.begin(), n);
	return buf.begin();
}

void File::read_c_str(std::string& s) {
	s = getdelim('\0');
}

void File::read_c_str(std::pmr::string& s) {
	s = getdelim('\0');
}

const char* File::getline() {
	line_buf_[0] = 0;
	ssize_t n = decompressor_->getline(&line_buf_, &line_buf_size_, file_);
	if (ferror(file_))
		throw runtime_error("Error reading file " + file_name_ + ". " + strerror(errno));
	while (n > 0) {
		if (line_buf_[n - 1] == '\n' || line_buf_[n - 1] == '\r')
			line_buf_[--n] = 0;
		else
			break;
	}
	++line_count;
	return line_buf_;
}

const char* File::getdelim(char delimiter) {
	line_buf_[0] = 0;
	ssize_t n = decompressor_->getdelim(&line_buf_, &line_buf_size_, delimiter, file_);
	if (ferror(file_))
		throw runtime_error("Error reading file " + file_name_ + ". " + strerror(errno));
	if (n > 0 && line_buf_[n - 1] == delimiter)
		line_buf_[n - 1] = 0;
	++line_count;
	return line_buf_;
}

size_t File::read_raw(string& dst, size_t count) {
	const size_t offset = dst.size();
	dst.resize(offset + count);
	const size_t n = read_max(&dst[offset], count);
	dst.resize(offset + n);
	return n;
}

uint64_t File::read_to_fasta_record_end(std::string& dst) {
	uint64_t n = 0;
	int c;
	for (;;) {
		if ((c = decompressor_->fgetc(file_)) == EOF)
			return n;
		if (c == '>' && (dst.back() == '\n' || dst.back() == '\r')) {
			decompressor_->ungetc(c, file_);
			return n;
		}
		++n;
		dst.push_back(c);
	}
}

bool File::eof() const {
	return feof(file_) != 0;
}

string File::peek(int64_t n) {	
	if (seekable_ && (decompressor_ == nullptr || decompressor_->lib() == CompressionLib::NONE)) {
		char buf[4];
		const size_t l = read_max(buf, n);
		seek(-(int64_t)l, SEEK_CUR);
		return string(buf, l);
	}
	else {
		string s;
		for (int64_t i = 0; i < n; ++i) {
			int c = decompressor_ ? decompressor_->fgetc(file_) : fgetc(file_);
			if (c == EOF)
				break;
			s.push_back((char)c);
		}
		for (auto i = s.rbegin(); i != s.rend(); ++i) {
			const int u = decompressor_ ? decompressor_->ungetc(*i, file_) : ungetc(*i, file_);
			if (u == EOF)
				throw std::runtime_error("Error calling ungetc on stream " + file_name_);
		}
		return s;
	}
}

int64_t File::count_lines(const std::string& file_name) {
	File f(file_name, "rb", File::Flags::DETECT_COMPRESSION);
	int64_t n = 0;
	const char* line;
	while (line = f.getline(), line[0] != '\0' || !f.eof())
		++n;
	return n;
}