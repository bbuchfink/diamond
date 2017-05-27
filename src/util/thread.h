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

#ifndef THREAD_H_
#define THREAD_H_

#include <vector>
#include <exception>
#include <stdexcept>
#include "fast_mutex.h"
#include "tinythread.h"
#include "util.h"

using tthread::thread;
using std::vector;

template<typename _t>
struct Atomic
{
	Atomic(const _t &v):
		v_ (v)
	{ }
	Atomic& operator=(const _t &v)
	{
		v_ = v;
		return *this;
	}
	volatile _t operator++(int)
	{
		mtx_.lock();
		_t r = v_++;
		mtx_.unlock();
		return r;
	}
	_t operator--(int)
	{
		mtx_.lock();
		_t r = v_--;
		mtx_.unlock();
		return r;
	}
private:
	volatile _t v_;
	tthread::mutex mtx_;
};

template<typename _context>
struct Thread_p
{
	Thread_p(unsigned thread_id, _context &context):
		thread_id (thread_id),
		context (&context)
	{ }
	unsigned thread_id;
	_context *context;
};

template<typename _context>
void pool_worker(void *p)
{
	((Thread_p<_context>*)p)->context->operator()(((Thread_p<_context>*)p)->thread_id);
	TLS::clear();
}

template<typename _context>
void launch_thread_pool(_context &context, unsigned threads)
{
	vector<tthread::thread*> t;
	vector<Thread_p<_context> > p;
	p.reserve(threads);
	for(unsigned i=0;i<threads;++i) {
		p.push_back(Thread_p<_context> (i, context));
		t.push_back(new tthread::thread(pool_worker<_context>, (void*)&p.back()));
	}
	for (vector<tthread::thread*>::iterator i = t.begin(); i != t.end(); ++i) {
		(*i)->join();
		delete *i;
	}
}

template<typename _context>
struct Schedule_context
{
	Schedule_context(_context &context, unsigned count):
		context (context),
		n (0),
		count (count)
	{ }
	void operator()(unsigned thread_id)
	{
		unsigned idx;
		while((idx = n++) < count)
			context(thread_id, idx);
	}
	_context &context;
	Atomic<unsigned> n;
	const unsigned count;
};

template<typename _context>
void launch_scheduled_thread_pool(_context &context, unsigned count, unsigned threads)
{
	Schedule_context<_context> c (context, count);
	launch_thread_pool(c, threads);
}

template<typename _f>
struct Thread_p0
{
	Thread_p0(_f f) :
		f(f)
	{ }
	_f f;
};

template<typename _f>
void thread_worker(void *p)
{
	Thread_p0<_f> *q = (Thread_p0<_f>*)p;
	q->f();
	delete q;
	TLS::clear();
}

template<typename _f>
thread* launch_thread(_f f)
{
	return new thread(thread_worker<_f>, new Thread_p0<_f>(f));
}


template<typename _f, typename _t1>
struct Thread_p1
{
	Thread_p1(_f f, _t1 p1) :
		f(f),
		p1(p1)
	{ }
	_f f;
	_t1 p1;
};

template<typename _f, typename _t1>
void thread_worker(void *p)
{
	Thread_p1<_f, _t1> *q = (Thread_p1<_f, _t1>*)p;
	q->f(q->p1);
	delete q;
	TLS::clear();
}

template<typename _f, typename _t1>
thread* launch_thread(_f f, _t1 p1)
{
	return new thread(thread_worker<_f, _t1>, new Thread_p1<_f, _t1>(f, p1));
}

template<typename _f, typename _t1, typename _t2>
struct Thread_p2
{
	Thread_p2(_f f, _t1 p1, _t2 p2) :
		f(f),
		p1(p1),
		p2(p2)
	{ }
	_f f;
	_t1 p1;
	_t2 p2;
};

template<typename _f, typename _t1, typename _t2>
void thread_worker(void *p)
{
	Thread_p2<_f, _t1, _t2> *q = (Thread_p2<_f, _t1, _t2>*)p;
	q->f(q->p1, q->p2);
	delete q;
	TLS::clear();
}

template<typename _f, typename _t1, typename _t2>
thread* launch_thread(_f f, _t1 p1, _t2 p2)
{
	return new thread(thread_worker<_f, _t1, _t2>, new Thread_p2<_f, _t1, _t2>(f, p1, p2));
}

template<typename _f, typename _t1, typename _t2, typename _t3>
struct Thread_p3
{
	Thread_p3(_f f, _t1 p1, _t2 p2, _t3 p3) :
		f(f),
		p1(p1),
		p2(p2),
		p3(p3)
	{ }
	_f f;
	_t1 p1;
	_t2 p2;
	_t3 p3;
};

template<typename _f, typename _t1, typename _t2, typename _t3>
void thread_worker(void *p)
{
	Thread_p3<_f, _t1, _t2, _t3> *q = (Thread_p3<_f, _t1, _t2, _t3>*)p;
	q->f(q->p1, q->p2, q->p3);
	delete q;
	TLS::clear();
}

template<typename _f, typename _t1, typename _t2, typename _t3>
thread* launch_thread(_f f, _t1 p1, _t2 p2, _t3 p3)
{
	return new thread(thread_worker<_f, _t1, _t2, _t3>, new Thread_p3<_f, _t1, _t2, _t3>(f, p1, p2, p3));
}

template<typename _f, typename _t1, typename _t2, typename _t3, typename _t4>
struct Thread_p4
{
	Thread_p4(_f f, _t1 p1, _t2 p2, _t3 p3, _t4 p4):
		f (f),
		p1 (p1),
		p2 (p2),
		p3 (p3),
		p4 (p4)
	{ }
	_f f;
	_t1 p1;
	_t2 p2;
	_t3 p3;
	_t4 p4;
};

template<typename _f, typename _t1, typename _t2, typename _t3, typename _t4>
void thread_worker(void *p)
{
	Thread_p4<_f,_t1,_t2,_t3,_t4> *q = (Thread_p4<_f,_t1,_t2,_t3,_t4>*)p;
	q->f(q->p1, q->p2, q->p3, q->p4);
	delete q;
	TLS::clear();
}

template<typename _f, typename _t1, typename _t2, typename _t3, typename _t4>
thread* launch_thread(_f f, _t1 p1, _t2 p2, _t3 p3, _t4 p4)
{ return new thread (thread_worker<_f,_t1,_t2,_t3,_t4>, new Thread_p4<_f,_t1,_t2,_t3,_t4> (f, p1, p2, p3, p4)); }

template<typename _f, typename _t1, typename _t2, typename _t3, typename _t4, typename _t5>
struct Thread_p5
{
	Thread_p5(_f f, _t1 p1, _t2 p2, _t3 p3, _t4 p4, _t5 p5) :
		f(f),
		p1(p1),
		p2(p2),
		p3(p3),
		p4(p4),
		p5(p5)
	{ }
	_f f;
	_t1 p1;
	_t2 p2;
	_t3 p3;
	_t4 p4;
	_t5 p5;
};

template<typename _f, typename _t1, typename _t2, typename _t3, typename _t4, typename _t5>
void thread_worker(void *p)
{
	Thread_p5<_f, _t1, _t2, _t3, _t4, _t5> *q = (Thread_p5<_f, _t1, _t2, _t3, _t4, _t5>*)p;
	q->f(q->p1, q->p2, q->p3, q->p4, q->p5);
	delete q;
	TLS::clear();
}

template<typename _f, typename _t1, typename _t2, typename _t3, typename _t4, typename _t5>
thread* launch_thread(_f f, _t1 p1, _t2 p2, _t3 p3, _t4 p4, _t5 p5)
{
	return new thread(thread_worker<_f, _t1, _t2, _t3, _t4, _t5>, new Thread_p5<_f, _t1, _t2, _t3, _t4, _t5>(f, p1, p2, p3, p4, p5));
}

template<typename _f, typename _t1, typename _t2, typename _t3, typename _t4, typename _t5, typename _t6>
struct Thread_p6
{
	Thread_p6(_f f, _t1 p1, _t2 p2, _t3 p3, _t4 p4, _t5 p5, _t6 p6) :
		f(f),
		p1(p1),
		p2(p2),
		p3(p3),
		p4(p4),
		p5(p5),
		p6(p6)
	{ }
	_f f;
	_t1 p1;
	_t2 p2;
	_t3 p3;
	_t4 p4;
	_t5 p5;
	_t6 p6;
};

template<typename _f, typename _t1, typename _t2, typename _t3, typename _t4, typename _t5, typename _t6>
void thread_worker(void *p)
{
	Thread_p6<_f, _t1, _t2, _t3, _t4, _t5,_t6> *q = (Thread_p6<_f, _t1, _t2, _t3, _t4, _t5,_t6>*)p;
	q->f(q->p1, q->p2, q->p3, q->p4, q->p5,q->p6);
	delete q;
	TLS::clear();
}

template<typename _f, typename _t1, typename _t2, typename _t3, typename _t4, typename _t5,typename _t6>
thread* launch_thread(_f f, _t1 p1, _t2 p2, _t3 p3, _t4 p4, _t5 p5, _t6 p6)
{
	return new thread(thread_worker<_f, _t1, _t2, _t3, _t4, _t5, _t6>, new Thread_p6<_f, _t1, _t2, _t3, _t4, _t5, _t6>(f, p1, p2, p3, p4, p5, p6));
}

struct Thread_pool : public vector<thread*>
{
	void join_all()
	{
		for (iterator i = begin(); i != end(); ++i) {
			(*i)->join();
			delete *i;
		}
		clear();
	}
};

#endif /* THREAD_H_ */
