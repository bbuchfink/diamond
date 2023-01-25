#include <queue>
#include "../io/serialize.h"
#include "../io/exceptions.h"

template<typename T, typename It, typename F>
void merge_sorted_files(const It begin, const It end, F& f) {
	struct Entry {
		Entry(T&& value, ptrdiff_t idx):
			value(value),
			idx(idx)
		{}
		Entry(const Entry& e):
			value(e.value),
			idx(e.idx)
		{}
		Entry& operator=(const Entry& e) {
			value = e.value;
			idx = e.idx;
			return *this;
		}
		bool operator<(const Entry& e) const {
			return !(value < e.value);
		}
		T value;
		ptrdiff_t idx;
	};
	std::priority_queue<Entry> q;
	std::vector<TypeDeserializer<T>> d;
	for (It i = begin; i != end; ++i) {
		d.emplace_back(*i);
		try {
			q.emplace(d.back().get(), i - begin);
		}
		catch (EndOfStream&) {
		}
	}
	while (!q.empty()) {
		f(q.top().value);
		const ptrdiff_t idx = q.top().idx;
		q.pop();
		try {
			q.emplace(d[idx].get(), idx);
		}
		catch (EndOfStream&) {
		}
	}
}