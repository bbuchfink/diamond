/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#include <chrono>
#include <thread>
#include "output.h"
#include "util/util.h"
#include "util/parallel/thread_pool.h"
#include "data/block/block.h"

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