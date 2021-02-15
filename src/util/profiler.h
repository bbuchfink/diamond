#pragma once
#include <map>
#include "log_stream.h"

struct Profiler {

	Profiler(const char* key) :
		timer(),
		key(key)
	{}

	~Profiler() {
		finish();
	}

	void finish() {
		if (!key) return;
		times[key] += timer.nanoseconds();
		key = nullptr;
	}

	static void print(size_t n) {
		for (const auto& i : times)
			message_stream << i.first << ": " << (double)i.second/n/1e3 << " micros" << std::endl;
	}

	task_timer timer;
	const char* key;
	static std::map<std::string, uint64_t> times;

};