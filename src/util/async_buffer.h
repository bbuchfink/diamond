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

#ifndef ASYNC_BUFFER_H_
#define ASYNC_BUFFER_H_

#include <vector>
#include <boost/thread.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

namespace io = boost::iostreams;

using std::vector;
using std::string;
using std::endl;
using boost::lockfree::queue;
using boost::thread;
using boost::ptr_vector;

struct Buffer_file_read_exception : public diamond_exception
{
	Buffer_file_read_exception(const char* file_name, size_t count, size_t n):
		diamond_exception (string("Error reading buffer file ") + file_name + " (" + boost::lexical_cast<string>(count) + '/' + boost::lexical_cast<string>(n) + ')')
	{ }
};

const unsigned async_buffer_max_bins = 4;

template<typename _t>
struct Async_buffer
{

	typedef vector<_t> Vector;

	Async_buffer(size_t input_count, const string &tmpdir, unsigned bins):
		bins_ (bins),
		bin_size_ ((input_count + bins_ - 1) / bins_),
		done_ (false),
		push_count_ (0)
	{
		log_stream << "Async_buffer() " << input_count << ',' << bin_size_ << endl;
		for(unsigned i=0;i<bins;++i) {
			tmp_file_.push_back(Temp_file ());
			out_.push_back(new Output_stream (tmp_file_[i]));
		}
		memset(size_, 0, sizeof(size_));
		writer_thread_  = new thread(writer, this);
	}

	struct Iterator
	{
		Iterator(Async_buffer &parent):
			parent_ (&parent)
		{
			for(unsigned i=0;i<parent_->bins_;++i)
				buffer_[i] = new vector<_t>;
		}
		void push(const _t &x)
		{
			//++parent_->push_count_;
			const unsigned bin = x / parent_->bin_size_;
			assert(bin < parent_->bins());
			buffer_[bin]->push_back(x);
			if(buffer_[bin]->size() == buffer_size) {
				parent_->out_queue_[bin].push(buffer_[bin]);
				buffer_[bin] = new vector<_t>;
			}
		}
		~Iterator()
		{
			for(unsigned bin=0;bin<parent_->bins_;++bin) {
				//log_stream << buffer_[bin]->size() << endl;
				//parent_->push_count_ += buffer_[bin]->size();
				parent_->out_queue_[bin].push(buffer_[bin]);
			}
		}
	private:
		vector<_t> *buffer_[async_buffer_max_bins];
		Async_buffer *parent_;
	};

	void close()
	{
		done_ = true;
		writer_thread_->join();
		delete writer_thread_;
		out_.clear();
		for(unsigned i=0;i<bins_;++i)
			log_stream << "Queue " << i << " status " << out_queue_[i].empty() << endl;
		log_stream << "Async_buffer.close() " << push_count_ << endl;
	}

	void load(vector<_t> &data, unsigned i) const
	{
		log_stream << "Async_buffer.load() " << size_[i] << endl;
		data.resize(size_[i]);
		if(size_[i] > 0) {
			Input_stream f (tmp_file_[i]);
			const size_t n = f.read(data.data(), size_[i]);
			f.close();
			if(n != size_[i])
				throw Buffer_file_read_exception(f.file_name.c_str(), size_[i], n);
		}
	}

	unsigned bins() const
	{ return bins_; }

private:

	enum { buffer_size = 65536 };

	void flush_queues()
	{
		vector<_t> *v;
		for(unsigned i=0;i<bins_;++i)
			while(out_queue_[i].pop(v)) {
				out_[i].write(v->data(), v->size());
				/*for(unsigned j=0;j<v->size();++j)
					cout << v->operator[](j);*/
				size_[i] += v->size();
				delete v;
			}
	}

	static void writer(Async_buffer *me)
	{
		try {
			while(!me->done_ && !exception_state()) {
				me->flush_queues();
				boost::this_thread::sleep_for(boost::chrono::milliseconds(1));
			}
			me->flush_queues();
		} catch(std::exception& e) {
			exception_state.set(e);
		}
	}

	const unsigned bins_, bin_size_;
	ptr_vector<Output_stream> out_;
#ifdef NDEBUG
	queue<vector<_t>*> out_queue_[async_buffer_max_bins];
#else
	queue<vector<_t>*,boost::lockfree::capacity<128> > out_queue_[async_buffer_max_bins];
#endif
	size_t size_[async_buffer_max_bins];
	thread *writer_thread_;
	bool done_;
	vector<Temp_file> tmp_file_;
	boost::atomic<size_t> push_count_;

};

#endif /* ASYNC_BUFFER_H_ */
