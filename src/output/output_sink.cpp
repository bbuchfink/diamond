/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

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

#include <chrono>
#include <thread>
#include "output.h"
#include "../data/queries.h"
#include "../util/util.h"
#include "../util/parallel/thread_pool.h"

using std::chrono::high_resolution_clock;
using std::chrono::seconds;
using std::chrono::milliseconds;
using std::chrono::duration_cast;
using std::endl;
using std::string;
using std::vector;

std::unique_ptr<ReorderQueue<TextBuffer*, OutputWriter>> output_sink;

void heartbeat_worker(size_t qend, const Search::Config* cfg)
{
	static const int interval = 100;
	static thread_local high_resolution_clock::time_point t0 = high_resolution_clock::now();
	int n = 0;
	size_t next;
	while ((next = output_sink->next()) < qend) {
		if (n == interval) {
			const string title(cfg->query->ids()[next]);
			verbose_stream << "Queries=" << next
				<< " size=" << megabytes(output_sink->size())
				<< " max_size=" << megabytes(output_sink->max_size())
				<< " next=" << title.substr(0, title.find(' '))
				<< " queue=" << cfg->thread_pool->queue_len(0) << "/" << cfg->thread_pool->queue_len(1)
				//<< " ETA=" << (double)duration_cast<seconds>(high_resolution_clock::now() - t0).count() / (next - OutputSink::get().begin()) * (qend - next) << "s"
				<< endl;
			n = 0;
		}
		else
			++n;
		std::this_thread::sleep_for(milliseconds(10));
	}
}