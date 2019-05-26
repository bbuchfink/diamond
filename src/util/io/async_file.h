#ifndef ASYNC_FILE_H_
#define ASYNC_FILE_H_

#include <mutex>
#include "temp_file.h"

struct AsyncFile : public TempFile {

	template<typename _t>
	void write(const _t *ptr, size_t count)
	{
		std::lock_guard<std::mutex> guard(mtx_);
		write_raw((const char*)ptr, count * sizeof(_t));
	}

private:

	std::mutex mtx_;

};

#endif