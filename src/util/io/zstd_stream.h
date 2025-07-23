#pragma once
#include <vector>
#include <zstd.h>
#include "stream_entity.h"
#include "input_stream_buffer.h"

struct ZstdSink : public StreamEntity
{
	ZstdSink(StreamEntity* prev);
	virtual void close() override;
	virtual void write(const char* ptr, size_t count) override;
private:
	ZSTD_CStream* stream;
};

struct ZstdSource : public StreamEntity
{
	ZstdSource(InputStreamBuffer* prev);
	virtual size_t read(char* ptr, size_t count) override;
	virtual void close() override;
	virtual void rewind() override;
	virtual bool eof() override;
private:
	void init();
	ZSTD_DStream* stream;
	ZSTD_inBuffer in_buf;
	bool eos_;
};
