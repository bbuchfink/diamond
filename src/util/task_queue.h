/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef TASK_QUEUE_H_
#define TASK_QUEUE_H_

#include "tinythread.h"

// #define ENABLE_LOGGING

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
