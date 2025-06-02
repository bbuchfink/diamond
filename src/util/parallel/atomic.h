#pragma once
#include "filestack.h"
#include "util/log_stream.h"

struct Atomic {

	Atomic(const std::string& file_name):
		stack_(file_name)
	{
	}

	int fetch_add(int n = 1) {
		std::string s;
		stack_.lock();
		int i;
		stack_.pop(i);
		if (i == -1)
			i = 0;
		stack_.push(i + n);
		// message_stream << "Atomic = " << i << std::endl;
		stack_.unlock();
		return i;
	}

	int get() {
		int i;
		stack_.top(i);
		return i >= 0 ? i : 0;
	}

	void await(int n) {
		for (;;) {
			int i;
			stack_.top(i);
			if (i >= n)
				return;
			std::this_thread::sleep_for(std::chrono::seconds(1));
		}
	}

private:
	FileStack stack_;

};