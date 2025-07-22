#pragma once
#ifdef WITH_ZSTD
#include <zstd.h>
#endif

struct CompressedBuffer {
	CompressedBuffer() :
		buf_(32768),
#ifdef WITH_ZSTD
		stream_(ZSTD_createCStream()),
#endif
		size_(0)
	{
#ifdef WITH_ZSTD
		if (!stream_)
			throw std::runtime_error("ZSTD_createCStream error");
#else
		throw std::runtime_error("ZSTD support missing. Compile with ZSTD library (CMake option -DWITH_ZSTD=ON).");
#endif
	}
	~CompressedBuffer() {
#ifdef WITH_ZSTD
		ZSTD_freeCStream(stream_);
#endif
	}
	void write(const char* ptr, int64_t n) {
#ifdef WITH_ZSTD
		ZSTD_inBuffer in_buf;
		in_buf.src = ptr;
		in_buf.size = n;
		in_buf.pos = 0;
		ZSTD_outBuffer out_buf;
		do {
			out_buf.dst = buf_.data();
			out_buf.size = buf_.size();
			out_buf.pos = size_;
			if (ZSTD_isError(ZSTD_compressStream(stream_, &out_buf, &in_buf)))
				throw std::runtime_error("ZSTD_compressStream");
			size_ = out_buf.pos;
			if (in_buf.pos < in_buf.size)
				buf_.resize(buf_.size() + 32768);
		} while (in_buf.pos < in_buf.size);
#endif
	}
	void finish() {
#ifdef WITH_ZSTD
		ZSTD_outBuffer out_buf;
		size_t n;
		do {
			out_buf.dst = buf_.data();
			out_buf.size = buf_.size();
			out_buf.pos = size_;
			if (ZSTD_isError(n = ZSTD_endStream(stream_, &out_buf)))
				throw std::runtime_error("ZSTD_endStream");
			size_ = out_buf.pos;
			if (n > 0)
				buf_.resize(buf_.size() + 32768);
		} while (n > 0);
		ZSTD_freeCStream(stream_);
		stream_ = nullptr;
#endif
	}
	void clear() {
#ifdef WITH_ZSTD
		stream_ = ZSTD_createCStream();
		size_ = 0;
#endif
	}
	const char* data() const {
		return buf_.data();
	}
	int64_t size() const {
		return size_;
	}
private:
	std::vector<char> buf_;
#ifdef WITH_ZSTD
	ZSTD_CStream* stream_;
#endif
	int64_t size_;
};