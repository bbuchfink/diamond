#include "../../basic/sequence.h"
#include "../data_structures/bit_vector.h"

struct KmerFilter {

	KmerFilter(Sequence seq, int k);
	std::pair<double, double> covered(Sequence seq) const;

private:

	KmerFilter(int k, Loc count, BitVector&& table);

	const int k_;
	const Loc count_;
	const BitVector table_;

	template<int K> friend KmerFilter build(Sequence);

};