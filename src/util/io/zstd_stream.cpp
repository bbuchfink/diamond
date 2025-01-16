#include "zstd_stream.h"

using std::pair;
using std::runtime_error;

ZstdSink::ZstdSink(StreamEntity* prev) :
	StreamEntity(prev),
	stream(ZSTD_createCStream())
{
	if (!stream)
		throw runtime_error("ZSTD_createCStream error");
}

void ZstdSink::write(const char* ptr, size_t count)
{
	ZSTD_inBuffer in_buf;
	in_buf.src = ptr;
	in_buf.size = count;
	in_buf.pos = 0;
	ZSTD_outBuffer out_buf;
	do {
		pair<char*, char*> out = prev_->write_buffer();
		out_buf.dst = out.first;
		out_buf.size = out.second - out.first;
		out_buf.pos = 0;
		if(ZSTD_isError(ZSTD_compressStream(stream, &out_buf, &in_buf)))
			throw runtime_error("ZSTD_compressStream");
		prev_->flush(out_buf.pos);
	} while (in_buf.pos < in_buf.size);
}

void ZstdSink::close()
{
	if (!stream)
		return;
	ZSTD_outBuffer out_buf;
	size_t n;
	do {
		pair<char*, char*> out = prev_->write_buffer();
		out_buf.dst = out.first;
		out_buf.size = out.second - out.first;
		out_buf.pos = 0;
		if (ZSTD_isError(n = ZSTD_endStream(stream, &out_buf)))
			throw runtime_error("ZSTD_endStream");
		prev_->flush(out_buf.pos);
	} while (n > 0);
	ZSTD_freeCStream(stream);
	stream = nullptr;
	prev_->close();
}

ZstdSource::ZstdSource(InputStreamBuffer* prev):
StreamEntity(prev)
{
	init();
}


void ZstdSource::init()
{
	eos_ = false;
	stream = ZSTD_createDStream();
	in_buf.pos = 0;
	in_buf.size = 0;
	in_buf.src = nullptr;
}

size_t ZstdSource::read(char* ptr, size_t count)
{
	ZSTD_outBuffer out_buf;
	out_buf.dst = ptr;
	out_buf.pos = 0;
	out_buf.size = count;
	InputStreamBuffer* buf = static_cast<InputStreamBuffer*>(prev_);
	while (out_buf.pos < out_buf.size && !eos_) {
		if (in_buf.pos >= in_buf.size) {
			if (buf->begin == buf->end)
				buf->fetch();
			in_buf.pos = 0;
			in_buf.size = buf->end - buf->begin;
			in_buf.src = buf->begin;
			buf->begin += in_buf.size;
			if (in_buf.size == 0) {
				eos_ = true;
				break;
			}
		}
		if (ZSTD_isError(ZSTD_decompressStream(stream, &out_buf, &in_buf)))
			throw runtime_error("ZSTD_decompressStream");
	}
	return out_buf.pos;
}

bool ZstdSource::eof() {
	return eos_;
}

void ZstdSource::close()
{
	if (!stream)
		return;
	ZSTD_freeDStream(stream);
	stream = nullptr;
	prev_->close();
}

void ZstdSource::rewind()
{
	prev_->rewind();
	init();
}