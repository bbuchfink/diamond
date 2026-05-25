/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/
// SPDX-License-Identifier: GPL-3.0-or-later

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