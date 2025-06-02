#pragma once
#include <vector>
#include <zstd.h>
#include "stream_entity.h"
#include "input_stream_buffer.h"

struct ZstdSink : public StreamEntity
{
	ZstdSink(StreamEntity* prev);
	virtual void close();
	virtual void write(const char* ptr, size_t count);
private:
	ZSTD_CStream* stream;
};

struct ZstdSource : public StreamEntity
{
	ZstdSource(InputStreamBuffer* prev);
	virtual size_t read(char* ptr, size_t count);
	virtual void close();
	virtual void rewind();
	virtual bool eof() override;
private:
	void init();
	ZSTD_DStream* stream;
	ZSTD_inBuffer in_buf;
	bool eos_;
};

struct CompressedBuffer {
	CompressedBuffer();
	~CompressedBuffer();
	void write(const char* ptr, int64_t n);
	void finish();
	void clear();
	const char* data() const {
		return buf_.data();
	}
	int64_t size() const {
		return size_;
	}
private:
	std::vector<char> buf_;
	ZSTD_CStream* stream_;
	int64_t size_;
};