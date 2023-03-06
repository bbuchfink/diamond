#pragma once
#include <vector>
#include "table.h"

using std::vector;

namespace Util { namespace Tsv {

struct BuildHelper {
	BuildHelper(int64_t size) {
		static const double TOLERANCE = 1.1;
		buffer.reserve(size_t(size * TOLERANCE));
		limits.push_back(0);
	}
	void add(Table* t) {
		buffer.insert(buffer.end(), t->data_.begin(), t->data_.end());
		limits.insert(limits.end(), t->limits_.begin() + 1, t->limits_.end());
		counts.push_back((int32_t)t->size());
	}
	Table get(const Schema& schema) {
		const int64_t n = counts.size();
		if (n == 0)
			return Table(schema);
		vector<int64_t> offsets;
		offsets.reserve(n);
		offsets.push_back(0);
		vector<int64_t>::const_iterator it = limits.begin() + counts[0];
		for (int64_t i = 1; i < n; ++i) {
			offsets.push_back(*it + offsets.back());
			it += counts[i];
		}
		vector<int64_t>::iterator begin = limits.begin() + counts[0] + 1;
		for (int64_t i = 1; i < n; ++i) {
			vector<int64_t>::iterator end = begin + counts[i];
			const auto d = offsets[i];
			for (auto it = begin; it < end; ++it)
				*it += d;
			begin += counts[i];
		}
		return Table(schema, move(buffer), move(limits));
	}
	vector<char> buffer;
	vector<int64_t> limits;
	vector<int32_t> counts;
};

}}