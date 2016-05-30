/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef TASK_QUEUE_H_
#define TASK_QUEUE_H_

#include "tinythread.h"

//#define ENABLE_LOGGING

template<typename _t, typename _callback>
struct Task_queue
{

	Task_queue(size_t limit, _callback &callback):
		queue_ (limit),
		state_ (limit),
		head_ (0),
		tail_ (0),
		limit_ (limit),
		head_idx_ (0),
		queued_ (0),
		queued_size_ (0),
		at_end_ (false),
		callback_ (callback)
	{ }

	volatile bool waiting() const
	{ return tail_ - head_ >= limit_; }

	template<typename _init>
	volatile bool get(size_t &n, _t*& res, _init &init)
	{
		{
			mtx_.lock();
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue get() thread=" << tthread::thread::get_current_thread_id() << " waiting=" << waiting() << " head=" << head_ << " tail=" << tail_ << endl;
#endif
			while(waiting() && !at_end_)
				cond_.wait(mtx_);
			if(at_end_) {
#ifdef ENABLE_LOGGING
				log_stream << "Task_queue get() thread=" << tthread::thread::get_current_thread_id() << " quit" << endl;
#endif
				mtx_.unlock();
				return false;
			}
			n = tail_++;
			res = &slot(n);
			if(!init())
				at_end_ = true;
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue get() thread=" << tthread::thread::get_current_thread_id() << " n=" << n << endl;
#endif
			mtx_.unlock();
		}
		if(at_end_)
			cond_.notify_all();
		return true;
	}

	void wake_all()
	{ cond_.notify_all(); }

	volatile void push(size_t n)
	{
		mtx_.lock();
		if(n == head_) {
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue flush() thread=" << tthread::thread::get_current_thread_id() << " n=" << n << endl;
#endif
			mtx_.unlock();
			flush();
		} else {
#ifdef ENABLE_LOGGING
			message_stream << "Task_queue push() thread=" << tthread::thread::get_current_thread_id() << " n=" << n << " head=" << head_ << " size=" << queued_size_ << endl;
#endif
			state_[idx(n)] = true;
			++queued_;
			queued_size_ += slot(n).size();
			mtx_.unlock();
		}
	}

	volatile unsigned flush()
	{
		bool next = false;
		unsigned n = 0;
		do {
			callback_(slot(head_));
			{
				mtx_.lock();
#ifdef ENABLE_LOGGING
				log_stream << "Task_queue flush() thread=" << tthread::thread::get_current_thread_id() << " head=" << head_ << " waiting=" << tail_-head_ << "/" << limit_ << " size=" << queued_size_ << endl;
#endif
				state_[idx(head_)] = false;
				++head_;
				head_idx_ = (head_idx_ + 1) % limit_;
				if (state_[idx(head_)]) {
					next = true;
#ifdef ENABLE_LOGGING
					log_stream << "Task_queue n=" << queued_ << endl;
#endif
					--queued_;
					queued_size_ -= slot(head_).size();
				}
				else
					next = false;
				mtx_.unlock();
			}
			cond_.notify_one();
			++n;
		} while (next);
		return n;
	}

private:

	size_t idx(size_t n) const
	{ return (head_idx_ + n - head_)%limit_; }

	_t& slot(size_t n)
	{ return queue_[idx(n)]; }

	vector<_t> queue_;
	vector<bool> state_;
	tthread::mutex mtx_;
	tthread::condition_variable cond_;
	volatile size_t head_, tail_, limit_, head_idx_, queued_, queued_size_;
	bool at_end_;
	_callback &callback_;

};

#endif
