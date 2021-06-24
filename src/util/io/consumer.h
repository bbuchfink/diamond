#pragma once

struct Consumer {
	virtual void consume(const char *ptr, size_t n) = 0;
	virtual void finalize() {}
	virtual ~Consumer() = default;
};