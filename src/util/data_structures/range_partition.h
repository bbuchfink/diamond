#include <algorithm>
#include <utility>
#include <assert.h>
#include <limits.h>

template<int _n, typename _t>
struct RangePartition {

	RangePartition(const int *begin, int count, int end) :
		count_(1)
	{
		assert(count > 0);
		std::pair<int, int> b[_n];
		for (int i = 0; i < count; ++i)
			b[i] = { begin[i], i };
		std::sort(b, b + count);
		std::fill(mask_[0], mask_[0] + _n, std::numeric_limits<_t>::min());
		mask_[0][b[0].second] = 0;
		begin_[0] = b[0].first;
		for (int i = 1; i < count; ++i) {
			if (begin_[count_ - 1] < b[i].first) {
				begin_[count_] = b[i].first;
				std::copy(mask_[count_ - 1], mask_[count_ - 1] + _n, mask_[count_]);
				mask_[count_][b[i].second] = 0;
				++count_;
			}
			else
				mask_[count_ - 1][b[i].second] = 0;
		}
		begin_[count_] = end;
	}

	int begin(int i) const {
		return begin_[i];
	}

	int end(int i) const {
		return begin_[i + 1];
	}

	int count() const {
		return count_;
	}

	const _t* mask(int i) const {
		return mask_[i];
	}

private:

	int begin_[_n + 1];
	_t mask_[_n][_n];
	int count_;
	
};