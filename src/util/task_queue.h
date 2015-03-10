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

struct Void_callback
{
	void operator()() const
	{ }
};

template<typename _t, typename _callback>
struct Task_queue
{

	Task_queue(unsigned n, unsigned limit, _callback &callback):
		queue_ (limit),
		head_ (0),
		tail_ (0),
		limit_ (limit),
		end_ (n),
		head_idx_ (0),
		callback_ (callback)
	{
		log_stream << "Task_queue() n=" << end_ << " limit=" << limit_ << endl;
	}

	volatile bool waiting() const
	{ return tail_ - head_ >= limit_; }

	template<typename _init>
	volatile bool get(size_t &n, _init &init)
	{
		{
			boost::unique_lock<boost::mutex> lock (mtx_);
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue get() thread=" << omp_get_thread_num() << " waiting=" << waiting() << " head=" << head_ << " tail=" << tail_ << endl;
#endif
			while(waiting() && tail_ < end_ && !exception_state())
				cond_.wait(lock);
			if(exception_state() || tail_ >= end_) {
#ifdef ENABLE_LOGGING
				log_stream << "Task_queue get() thread=" << omp_get_thread_num() << " quit" << endl;
#endif
				return false;
			}
			n = tail_++;
			init();
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue get() thread=" << omp_get_thread_num() << " n=" << n << endl;
#endif
		}
		if(tail_ == end_)
			cond_.notify_all();
		return true;
	}

	void wake_all()
	{ cond_.notify_all(); }

	volatile void push(unsigned n, const _t& v)
	{
		mtx_.lock();
		if(n == head_) {
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue flush() thread=" << omp_get_thread_num() << " n=" << n << endl;
#endif
			mtx_.unlock();
			flush(v);
		} else {
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue push() thread=" << omp_get_thread_num() << "n=" << n << " head=" << head_ << endl;
#endif
			slot(n) = v;
			mtx_.unlock();
		}
	}

	volatile unsigned flush(const _t &v)
	{
		callback_(v);
		bool notify, next=false;
		_t nextv;
		{
			boost::lock_guard<boost::mutex> lock(mtx_);
			notify = waiting();
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue flush() thread=" << omp_get_thread_num() << " head=" << head_ << " waiting=" << notify << endl;
#endif
			slot(head_) = _t ();
			++head_;
			head_idx_ = (head_idx_+1)%limit_;
			if(slot(head_)) {
				nextv = slot(head_);
				next = true;
			}
		}
		if(notify)
			cond_.notify_one();
		if(next)
			return flush(nextv)+1;
		return 1;
	}

private:

	_t& slot(unsigned n)
	{ return queue_[(head_idx_ + n - head_)%limit_]; }

	vector<_t> queue_;
	boost::mutex mtx_;
	boost::condition_variable cond_;
	volatile unsigned head_, tail_, limit_, end_, head_idx_;
	_callback &callback_;

};

template<typename _t, typename _callback>
struct Task_queue2
{

	Task_queue2(unsigned n, unsigned limit, _callback &callback):
		queue_ (limit),
		state_ (limit),
		head_ (0),
		tail_ (0),
		limit_ (limit),
		end_ (n),
		head_idx_ (0),
		callback_ (callback)
	{
		log_stream << "Task_queue() n=" << end_ << " limit=" << limit_ << endl;
	}

	volatile bool waiting() const
	{ return tail_ - head_ >= limit_; }

	template<typename _init>
	volatile bool get(size_t &n, _t*& res, _init &init)
	{
		{
			boost::unique_lock<boost::mutex> lock (mtx_);
			log_stream << "Task_queue get() thread=" << omp_get_thread_num() << " waiting=" << waiting() << " head=" << head_ << " tail=" << tail_ << endl;
			while(waiting() && tail_ < end_ && !exception_state())
				cond_.wait(lock);
			if(exception_state() || tail_ >= end_) {
				log_stream << "Task_queue get() thread=" << omp_get_thread_num() << " quit" << endl;
				return false;
			}
			n = tail_++;
			res = &slot(n);
			init();
			log_stream << "Task_queue get() thread=" << omp_get_thread_num() << " n=" << n << endl;
		}
		if(tail_ == end_)
			cond_.notify_all();
		return true;
	}

	void wake_all()
	{ cond_.notify_all(); }

	volatile void push(unsigned n)
	{
		mtx_.lock();
		if(n == head_) {
			log_stream << "Task_queue flush() thread=" << omp_get_thread_num() << " n=" << n << endl;
			mtx_.unlock();
			//printf("%u\n",flush(v));
			flush();
		} else {
			log_stream << "Task_queue push() thread=" << omp_get_thread_num() << "n=" << n << " head=" << head_ << endl;
			state_[idx(n)] = true;
			mtx_.unlock();
		}
	}

	volatile unsigned flush()
	{
		callback_(slot(head_));
		bool notify, next=false;
		{
			boost::lock_guard<boost::mutex> lock(mtx_);
			notify = waiting();
			log_stream << "Task_queue flush() thread=" << omp_get_thread_num() << " head=" << head_ << " waiting=" << notify << endl;
			state_[idx(head_)] = false;
			++head_;
			head_idx_ = (head_idx_+1)%limit_;
			if(state_[idx(head_)])
				next = true;
		}
		if(notify)
			cond_.notify_one();
		if(next)
			return flush()+1;
		return 1;
	}

private:

	unsigned idx(unsigned n) const
	{ return (head_idx_ + n - head_)%limit_; }

	_t& slot(unsigned n)
	{ return queue_[idx(n)]; }

	vector<_t> queue_;
	vector<bool> state_;
	boost::mutex mtx_;
	boost::condition_variable cond_;
	volatile unsigned head_, tail_, limit_, end_, head_idx_;
	_callback &callback_;

};

template<typename _t, typename _callback>
struct Task_queue3
{

	Task_queue3(unsigned limit, _callback &callback):
		queue_ (limit),
		head_ (0),
		tail_ (0),
		limit_ (limit),
		head_idx_ (0),
		at_end_ (false),
		callback_ (callback)
	{ }

	volatile bool waiting() const
	{ return tail_ - head_ >= limit_; }

	template<typename _init>
	volatile bool get(size_t &n, _init &init)
	{
		{
			boost::unique_lock<boost::mutex> lock (mtx_);
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue get() thread=" << omp_get_thread_num() << " waiting=" << waiting() << " head=" << head_ << " tail=" << tail_ << endl;
#endif
			while(waiting() && !at_end_ && !exception_state())
				cond_.wait(lock);
			if(exception_state() || at_end_) {
#ifdef ENABLE_LOGGING
				log_stream << "Task_queue get() thread=" << omp_get_thread_num() << " quit" << endl;
#endif
				return false;
			}
			n = tail_++;
			if(!init())
				at_end_ = true;
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue get() thread=" << omp_get_thread_num() << " n=" << n << endl;
#endif
		}
		if(at_end_)
			cond_.notify_all();
		return true;
	}

	void wake_all()
	{ cond_.notify_all(); }

	volatile void push(unsigned n, const _t& v)
	{
		mtx_.lock();
		if(n == head_) {
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue flush() thread=" << omp_get_thread_num() << " n=" << n << endl;
#endif
			mtx_.unlock();
			flush(v);
		} else {
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue push() thread=" << omp_get_thread_num() << "n=" << n << " head=" << head_ << endl;
#endif
			slot(n) = v;
			mtx_.unlock();
		}
	}

	volatile unsigned flush(const _t &v)
	{
		callback_(v);
		bool notify, next=false;
		_t nextv;
		{
			boost::lock_guard<boost::mutex> lock(mtx_);
			notify = waiting();
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue flush() thread=" << omp_get_thread_num() << " head=" << head_ << " waiting=" << notify << endl;
#endif
			slot(head_) = _t ();
			++head_;
			head_idx_ = (head_idx_+1)%limit_;
			if(slot(head_)) {
				nextv = slot(head_);
				next = true;
			}
		}
		if(notify)
			cond_.notify_one();
		if(next)
			return flush(nextv)+1;
		return 1;
	}

private:

	_t& slot(unsigned n)
	{ return queue_[(head_idx_ + n - head_)%limit_]; }

	vector<_t> queue_;
	boost::mutex mtx_;
	boost::condition_variable cond_;
	volatile unsigned head_, tail_, limit_, head_idx_;
	bool at_end_;
	_callback &callback_;

};


#endif
