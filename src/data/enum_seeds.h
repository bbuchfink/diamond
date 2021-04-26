#pragma once

#include "sequence_set.h"

template<typename _f, typename _filter>
void enum_seeds(const SequenceSet* seqs, _f* f, unsigned begin, unsigned end, std::pair<size_t, size_t> shape_range, const _filter* filter)
{
	vector<Letter> buf(seqs->max_len(begin, end));
	uint64_t key;
	for (unsigned i = begin; i < end; ++i) {
		const Sequence seq = (*seqs)[i];
		Reduction::reduce_seq(seq, buf);
		for (size_t shape_id = shape_range.first; shape_id < shape_range.second; ++shape_id) {
			const Shape& sh = shapes[shape_id];
			if (seq.length() < sh.length_) continue;
			Seed_iterator it(buf, sh);
			size_t j = 0;
			while (it.good()) {
				if (it.get(key, sh))
					if (filter->contains(key, shape_id))
						(*f)(key, seqs->position(i, j), shape_id);
				++j;
			}
		}
	}
	f->finish();
}

template<typename _f, uint64_t _b, typename _filter>
void enum_seeds_hashed(const SequenceSet* seqs, _f* f, unsigned begin, unsigned end, std::pair<size_t, size_t> shape_range, const _filter* filter)
{
	uint64_t key;
	for (unsigned i = begin; i < end; ++i) {
		const Sequence seq = (*seqs)[i];
		for (size_t shape_id = shape_range.first; shape_id < shape_range.second; ++shape_id) {
			const Shape& sh = shapes[shape_id];
			if (seq.length() < sh.length_) continue;
			const uint64_t shape_mask = sh.long_mask();
			//const __m128i shape_mask = sh.long_mask_sse_;
			Hashed_seed_iterator<_b> it(seq, sh);
			size_t j = 0;
			while (it.good()) {
				if (it.get(key, shape_mask)) {
					if (filter->contains(key, shape_id))
						(*f)(key, seqs->position(i, j), shape_id);
				}
				++j;
			}
		}
	}
	f->finish();
}

template<typename _f, typename _filter>
static void enum_seeds_worker(_f* f, const SequenceSet* seqs, unsigned begin, unsigned end, std::pair<size_t, size_t> shape_range, const _filter* filter, bool hashed)
{
	if (hashed) {
		const uint64_t b = Reduction::reduction.bit_size();
		switch (b) {
		case 4:
			enum_seeds_hashed<_f, 4, _filter>(seqs, f, begin, end, shape_range, filter);
			break;
		default:
			throw std::runtime_error("Unsupported reduction.");
		}
	}
	else
		enum_seeds<_f, _filter>(seqs, f, begin, end, shape_range, filter);
}

struct No_filter
{
	bool contains(uint64_t seed, uint64_t shape) const
	{
		return true;
	}
};

extern No_filter no_filter;

template <typename _f, typename _filter>
void enum_seeds(const SequenceSet* seqs, PtrVector<_f>& f, const std::vector<size_t>& p, size_t shape_begin, size_t shape_end, const _filter* filter, bool hashed)
{
	std::vector<std::thread> threads;
	for (unsigned i = 0; i < f.size(); ++i)
		threads.emplace_back(enum_seeds_worker<_f, _filter>, &f[i], seqs, (unsigned)p[i], (unsigned)p[i + 1], std::make_pair(shape_begin, shape_end), filter, hashed);
	for (auto& t : threads)
		t.join();
}