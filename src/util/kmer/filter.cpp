#include "filter.h"
#include "kmer.h"

using std::pair;
using std::runtime_error;
using std::move;

template<int K>
static KmerFilter build(Sequence seq) {
	KmerIterator<K> it(seq);
	BitVector v(power(20, K));
	Loc n = 0;
	while (it.good()) {
		assert(*it < v.size());
		v.set(*it);
		++it;
		++n;
	}
	return KmerFilter{ K, n, move(v) };
}

static KmerFilter build(Sequence seq, int k) {
	switch (k) {
	case 2:
		return build<2>(seq);
	case 3:
		return build<3>(seq);
	case 4:
		return build<4>(seq);
	case 5:
		return build<5>(seq);
	default:
		throw runtime_error("Unsupported kmer size");
	}
}

KmerFilter::KmerFilter(int k, Loc count, BitVector&& table):
	k_(k),
	count_(count),
	table_(move(table))
{}

KmerFilter::KmerFilter(Sequence seq, int k) :
	KmerFilter(build(seq, k))
{
}

template<int K>
static pair<Loc, Loc> covered(const BitVector& v, Sequence seq) {
	KmerIterator<K> it(seq);
	Loc n = 0, s = 0;
	while (it.good()) {
		if (v.get(*it))
			++s;
		++it;
		++n;
	}
	return { n,s };
}

static pair<Loc, Loc> covered(const BitVector& v, Sequence seq, int k) {
	switch (k) {
	case 2:
		return covered<2>(v, seq);
	case 3:
		return covered<3>(v, seq);
	case 4:
		return covered<4>(v, seq);
	case 5:
		return covered<5>(v, seq);
	default:
		throw runtime_error("Unsupported kmer size");
	}
}

pair<double, double> KmerFilter::covered(Sequence seq) const {
	const pair<Loc, Loc> n = ::covered(table_, seq, k_);
	return { (double)n.second / count_, (double)n.second / n.first };
}