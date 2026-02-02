/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include "filestack.h"
#include "util/log_stream.h"

struct Atomic {

	Atomic(const std::string& file_name):
		stack_(file_name)
	{
	}

	int64_t fetch_add(int64_t n = 1) {
		std::string s;
		stack_.lock();
		int64_t i;
		stack_.pop(i);
		if (i == -1)
			i = 0;
		stack_.push(i + n);
		// message_stream << "Atomic = " << i << std::endl;
		stack_.unlock();
		return i;
	}

	int64_t get() {
		int64_t i;
		stack_.top(i);
		return i >= 0 ? i : 0;
	}

	void await(int64_t n) {
		for (;;) {
			int64_t i;
			stack_.top(i);
			if (i >= n)
				return;
			std::this_thread::sleep_for(std::chrono::seconds(1));
		}
	}

private:
	FileStack stack_;

};