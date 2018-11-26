#ifndef UTIL_IO_CONSUMER_H_
#define UTIL_IO_CONSUMER_H_

struct Consumer {
	virtual void consume(const char *ptr, size_t n) = 0;
	virtual void finalize() {}
	virtual ~Consumer() = default;
};

#endif