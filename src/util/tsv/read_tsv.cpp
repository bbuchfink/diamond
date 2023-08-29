#include <queue>
#include <condition_variable>
#include <mutex>
#include <thread>
#include "tsv.h"
#include "../util/io/input_stream_buffer.h"
#include "file.h"

using std::mutex;
using std::condition_variable;
using std::runtime_error;
using std::string;
using std::function;
using std::queue;
using std::vector;
using std::search;
using std::copy;
using std::unique_lock;
using std::move;
using std::thread;
using std::back_inserter;

namespace Util { namespace Tsv {

static const int64_t READ_SIZE = 1 * (1 << 20);

void File::read(int64_t max_size, int threads, function<void(int64_t chunk, const char*, const char*)>& callback) {
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
			int64_t n = file_->read_raw(buf.data(), READ_SIZE);
			if (n == READ_SIZE) {
				file_->read_to(back_inserter(buf), '\n');
				n = buf.size();
			}
			total += n;
			if (n > 0) {
				vector<char> new_buf(buf.begin(), buf.begin() + n);
				{
					unique_lock<mutex> lock(mtx);
					read_cv.wait(lock, [&buffers, consumers] { return buffers.size() < consumers; });
					buffers.push(move(new_buf));
				}
			}
			if (n < READ_SIZE || total + READ_SIZE > max_size) {
				stop = true;
				consume_cv.notify_all();
				break;
			} else
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
					buf = move(buffers.front());
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

}}