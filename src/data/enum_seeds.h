#pragma once

#include "sequence_set.h"

enum class SeedEncoding { SPACED_FACTOR, HASHED, CONTIGUOUS };

template<typename _f, typename _filter>
void enum_seeds(SequenceSet* seqs, _f* f, unsigned begin, unsigned end, std::pair<size_t, size_t> shape_range, const _filter* filter, const std::vector<bool>* skip)
{
	vector<Letter> buf(seqs->max_len(begin, end));
	uint64_t key;
	for (unsigned i = begin; i < end; ++i) {
		if (skip && (*skip)[i / align_mode.query_contexts])
			continue;
		seqs->convert_to_std_alph(i);
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
						(*f)(key, seqs->position(i, j), i, shape_id);
				++j;
			}
		}
	}
	f->finish();
}

template<typename _f, uint64_t _b, typename _filter>
void enum_seeds_hashed(SequenceSet* seqs, _f* f, unsigned begin, unsigned end, std::pair<size_t, size_t> shape_range, const _filter* filter, const std::vector<bool>* skip)
{
	uint64_t key;
	for (unsigned i = begin; i < end; ++i) {
		if (skip && (*skip)[i / align_mode.query_contexts])
			continue;
		seqs->convert_to_std_alph(i);
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
						(*f)(key, seqs->position(i, j), i, shape_id);
				}
				++j;
			}
		}
	}
	f->finish();
}

template<typename _f, typename _it, typename _filter>
void enum_seeds_contiguous(SequenceSet* seqs, _f* f, unsigned begin, unsigned end, const _filter* filter, const std::vector<bool>* skip)
{
	uint64_t key;
	for (unsigned i = begin; i < end; ++i) {
		if (skip && (*skip)[i / align_mode.query_contexts])
			continue;
		seqs->convert_to_std_alph(i);
		const Sequence seq = (*seqs)[i];
		if (seq.length() < _it::length()) continue;
		_it it(seq);
		size_t j = 0;
		while (it.good()) {
			if (it.get(key))
				if (filter->contains(key, 0))
					if ((*f)(key, seqs->position(i, j), i, 0) == false)
						return;
			++j;
		}
	}
	f->finish();
}

template<typename _f, typename _filter, typename IteratorFilter>
static void enum_seeds_worker(_f* f, SequenceSet* seqs, unsigned begin, unsigned end, std::pair<size_t, size_t> shape_range, const _filter* filter, SeedEncoding code, const std::vector<bool>* skip)
{
	static const char* errmsg = "Unsupported contiguous seed.";
	if (code == SeedEncoding::CONTIGUOUS) {
		const uint64_t b = Reduction::reduction.bit_size(), l = shapes[shape_range.first].length_;
		switch (l) {
		case 7:
			switch (b) {
			case 4:
				enum_seeds_contiguous<_f, Contiguous_seed_iterator<7, 4, IteratorFilter>, _filter>(seqs, f, begin, end, filter, skip);
				break;
			default:
				throw std::runtime_error(errmsg);
			}
			break;
		case 6:
			switch (b) {
			case 4:
				enum_seeds_contiguous<_f, Contiguous_seed_iterator<6, 4, IteratorFilter>, _filter>(seqs, f, begin, end, filter, skip);
				break;
			default:
				throw std::runtime_error(errmsg);
			}
			break;
		case 5:
			switch (b) {
			case 4:
				enum_seeds_contiguous<_f, Contiguous_seed_iterator<5, 4, IteratorFilter>, _filter>(seqs, f, begin, end, filter, skip);
				break;
			default:
				throw std::runtime_error(errmsg);
			}
			break;
		default:
			throw std::runtime_error(errmsg);
		}
	}
	else if (code == SeedEncoding::HASHED) {
		const uint64_t b = Reduction::reduction.bit_size();
		switch (b) {
		case 4:
			enum_seeds_hashed<_f, 4, _filter>(seqs, f, begin, end, shape_range, filter, skip);
			break;
		default:
			throw std::runtime_error("Unsupported reduction.");
		}
	}
	else
		enum_seeds<_f, _filter>(seqs, f, begin, end, shape_range, filter, skip);
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
void enum_seeds(SequenceSet* seqs, PtrVector<_f>& f, const std::vector<size_t>& p, size_t shape_begin, size_t shape_end, const _filter* filter, SeedEncoding code, const std::vector<bool>* skip, bool filter_masked_seeds)
{
	std::vector<std::thread> threads;
	for (unsigned i = 0; i < f.size(); ++i)
		if (filter_masked_seeds)
			threads.emplace_back(enum_seeds_worker<_f, _filter, FilterMaskedSeeds>, &f[i], seqs, (unsigned)p[i], (unsigned)p[i + 1], std::make_pair(shape_begin, shape_end), filter, code, skip);
		else
			threads.emplace_back(enum_seeds_worker<_f, _filter, void>, &f[i], seqs, (unsigned)p[i], (unsigned)p[i + 1], std::make_pair(shape_begin, shape_end), filter, code, skip);
	for (auto& t : threads)
		t.join();
	seqs->alphabet() = Alphabet::STD;
}