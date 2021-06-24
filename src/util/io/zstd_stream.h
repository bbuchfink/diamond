#pragma once
#include <zstd.h>
#include "stream_entity.h"

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
	ZstdSource(StreamEntity* prev);
	virtual size_t read(char* ptr, size_t count);
	virtual void close();
	virtual void rewind();
private:
	void init();
	ZSTD_DStream* stream;
	ZSTD_inBuffer in_buf;
	bool eos_;
};