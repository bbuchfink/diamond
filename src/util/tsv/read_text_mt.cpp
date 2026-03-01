/****
Copyright (C) 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

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

#include <queue>
#include <condition_variable>
#include <mutex>
#include <thread>
#include "util/io/text_input_file.h"
#include "file.h"

using std::mutex;
using std::condition_variable;
using std::function;
using std::queue;
using std::vector;
using std::unique_lock;
using std::thread;
using std::back_inserter;
using std::lock_guard;

static const int64_t READ_SIZE = 1 * (1 << 20);

void read_text_mt(TextInputFile& file, int64_t max_size, int threads, function<void(int64_t chunk, const char*, const char*)>& callback) {
	queue<vector<char>> buffers;
	int64_t next_chunk = 0;
	bool stop = false;
	const size_t consumers = std::max(threads - 1, 1);
	mutex mtx;
	condition_variable consume_cv, read_cv;

	auto reader = [&]() {
		vector<char> buf(READ_SIZE);
		int64_t total = 0;
		for (;;) {
			int64_t n = file.read_raw(buf.data(), READ_SIZE);
			if (n == READ_SIZE) {
				file.read_to(back_inserter(buf), '\n');
				n = buf.size();
			}
			total += n;
			if (n > 0) {
				vector<char> new_buf(buf.begin(), buf.begin() + n);
				{
					unique_lock<mutex> lock(mtx);
					read_cv.wait(lock, [&buffers, consumers] { return buffers.size() < consumers; });
					buffers.push(std::move(new_buf));
				}
			}
			if (n < READ_SIZE || total + READ_SIZE > max_size) {
				{
					lock_guard<mutex> lock(mtx);
					stop = true;
				}
				consume_cv.notify_all();
				break;
			}
			else
				consume_cv.notify_one();
			buf.resize(READ_SIZE);
		}
		};

	auto consumer = [&] {
		vector<char> buf;
		int64_t chunk;
		for (;;) {
			{
				unique_lock<mutex> lock(mtx);
				consume_cv.wait(lock, [&buffers, &stop] { return stop || !buffers.empty(); });
				if (!buffers.empty()) {
					buf = std::move(buffers.front());
					buffers.pop();
					chunk = next_chunk++;
				}
				else
					return;
			}
			read_cv.notify_one();
			callback(chunk, buf.data(), buf.data() + buf.size());
		}
		};

	vector<thread> thread;
	thread.emplace_back(reader);
	for (size_t i = 0; i < consumers; ++i)
		thread.emplace_back(consumer);
	for (auto& t : thread)
		t.join();
}