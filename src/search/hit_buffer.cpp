#include "hit_buffer.h"
#include "util/log_stream.h"
#include "basic/config.h"

using std::vector;
using std::string;
using std::thread;
using std::copy;
using std::runtime_error;
using std::mutex;

namespace Search {

static int64_t load_hits(MMap& f, bool long_subject_offsets, int threads, Hit* out) {
	mutex mtx;
	const uint8_t* input_ptr = f.data(), * input_end = input_ptr + f.size();
	std::atomic_int64_t next(0);
	//SingleProducerQueue<std::vector<uint8_t>> queue(256, threads);
	auto worker = [&]() {
		for (;;) {
			uint32_t buf_size;
			const uint8_t* ptr;
			{
				std::lock_guard<mutex> lock(mtx);
				if(input_ptr >= input_end)
					return;
				copy(input_ptr, input_ptr + sizeof(buf_size), (char*)&buf_size);
				ptr = input_ptr + sizeof(buf_size);
				input_ptr += sizeof(buf_size) + buf_size;
			}
			//vector<uint8_t> buf;
			//buf.reserve(buf_size);
			//copy(buf_ptr, buf_ptr + buf_size, std::back_inserter(buf));
			std::vector<Hit> dst;
			//const uint8_t* ptr = buf.data(), * end = ptr + buf.size();
			const uint8_t* end = ptr + buf_size;
			uint16_t x;
			std::copy(ptr, ptr + sizeof(x), (char*)&x);
			ptr += sizeof(x);
			assert(x == 0);
			for (;;) {
				if (ptr >= end)
					break;
				uint32_t query_id, seed_offset;
				std::copy(ptr, ptr + sizeof(query_id), (char*)&query_id);
				ptr += sizeof(query_id);
				std::copy(ptr, ptr + sizeof(seed_offset), (char*)&seed_offset);
				ptr += sizeof(seed_offset);
				PackedLoc subject_loc;
				uint32_t x;
				for (;;) {
					if (ptr >= end)
						break;
					uint16_t score;
					std::copy(ptr, ptr + sizeof(score), (char*)&score);
					ptr += sizeof(score);
					if (score == 0)
						break;
					if (long_subject_offsets) {
						std::copy(ptr, ptr + sizeof(subject_loc), (char*)&subject_loc);
						ptr += sizeof(subject_loc);
					}
					else {
						std::copy(ptr, ptr + sizeof(uint32_t), (char*)&x);
						ptr += sizeof(uint32_t);
						subject_loc = x;
					}
#ifdef HIT_KEEP_TARGET_ID
					uint32_t target_block_id;
					std::copy(ptr, ptr + sizeof(target_block_id), (char*)&target_block_id);
					ptr += sizeof(target_block_id);
					dst.emplace_back(query_id, subject_loc, (Loc)seed_offset, score, target_block_id);
#else
					dst.emplace_back(query_id, subject_loc, (Loc)seed_offset, score);
#endif
				}
			}
			assert(!dst.empty());
			std::copy(dst.begin(), dst.end(), out + next.fetch_add(dst.size(), std::memory_order_relaxed));
		}
		};
	vector<thread> workers;
	for (int i = 0; i < threads; ++i)
		workers.emplace_back(worker);
	//int64_t max_queue_size = 0;
	//while (ptr < end) {
	//	queue.push(std::move(buf));
	//}
	//queue.close();
	for (auto& t : workers) {
		if (t.joinable())
			t.join();
	}
	//assert(queue.empty());
	return next;
}

HitBuffer::HitBuffer(const vector<Key>& key_partition, const string& tmpdir, bool long_subject_offsets, int query_contexts) :
	key_partition_(key_partition),
	long_subject_offsets_(long_subject_offsets),
	query_contexts_(query_contexts),
	bins_processed_(0),
	total_disk_size_(0),
	bin_mutex_(key_partition.size()),
	load_worker_(nullptr)
{
	log_stream << "Async_buffer() " << key_partition.back() << std::endl;
	count_ = new std::atomic_size_t[key_partition.size()];
	MMap::Options options;
	options.temp_dir_override = tmpdir;
	for (size_t i = 0; i < key_partition.size(); ++i) {
		if (config.trace_pt_membuf) {
			hit_buf_.emplace_back();
			hit_buf_.back().reserve(1ll * 1024ll * 1024ll * 1024ll / sizeof(Hit));
		}
		else
			tmp_file_.emplace_back(options);
		count_[i] = (size_t)0;
	}
}

void HitBuffer::load(int64_t max_size) {
	max_size = std::max(max_size, (int64_t)1);
	auto worker = [this](int end) {
		Hit* out = data_next_;
		for (; bins_processed_ < end; ++bins_processed_) {
			load_bin(out, bins_processed_);
			out += count_[bins_processed_];
		}
		};
	if (bins_processed_ == bins()) {
		data_next_ = nullptr;
		return;
	}
	int64_t size = count_[bins_processed_], current_size;
	int end = bins_processed_ + 1;
	if (!config.trace_pt_membuf) {
		int64_t disk_size = tmp_file_[bins_processed_].size();
		while (end < bins() && (size + (current_size = count_[end])) * ENTRY_SIZE < max_size && (end - bins_processed_ == 0)) {
			size += current_size;
			disk_size += tmp_file_[end].size();
			++end;
		}
		log_stream << "Async_buffer.load() " << size << "(" << (double)size * sizeof(Hit) / (1 << 30) << " GB, " << (double)disk_size / (1 << 30) << " GB on disk)" << std::endl;
		total_disk_size_ += disk_size;
		data_size_next_ = size;
		load_worker_ = new std::thread(worker, end);
	}
	input_range_next_.first = begin(bins_processed_);
	input_range_next_.second = this->end(end - 1);	
}

void HitBuffer::load_bin(Hit* out, int bin)
{
	if (!config.trace_pt_membuf) {
#ifndef _MSC_VER
		if (bin < bins() - 1)
			posix_madvise(tmp_file_[bin + 1].data(), tmp_file_[bin + 1].size(), POSIX_MADV_SEQUENTIAL | POSIX_MADV_WILLNEED);
#endif
		if (count_[bin] > 0) {
			//mprotect(out, )
			const size_t n = load_hits(tmp_file_[bin], long_subject_offsets_, 8, out);
			if (n != count_[bin])
				throw runtime_error("Mismatching hit count / possibly corrupted temporary file");
		}
		tmp_file_[bin].close();
	}
}

}