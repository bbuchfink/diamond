#ifndef EDGE_VEC_H_
#define EDGE_VEC_H_

#include <unordered_map>
#include <string>
#include <array>
#include <math.h>
#include <algorithm>
#include <vector>
#include "../io/temp_file.h"

namespace Util { namespace Algo { namespace UPGMA_MC {

enum class DistType { EVALUE, BITSCORE };

DistType dist_type();

struct CompactEdge {
	bool operator<(const CompactEdge &e) const {
		return d < e.d;
	}
	operator bool() const {
		return n1 != 0 || n2 != 0;
	}
	int n1, n2;
	double d;
};

struct EdgeVec {
	static constexpr int BUCKET_COUNT = 330;
	EdgeVec(const char *file);
	size_t nodes() const {
		return node_count_ > 0 ? node_count_ : acc2idx.size();
	}
	size_t size() const {
		return size_;
	}
	CompactEdge get();
	std::string print(int idx) const;
private:
	static int bucket(double d, DistType dt) {
		if (dt == DistType::EVALUE)
			return d == 0.0 ? 0 : std::min(323 + (int)std::log10(d), BUCKET_COUNT - 1);
		else
			return BUCKET_COUNT - 1 - std::max(std::min(int(-d * 300.0), BUCKET_COUNT - 1), 0);
	}	
	std::unordered_map<std::string, int> acc2idx;
	std::unordered_map<int, std::string> idx2acc;
	std::array<TempFile, BUCKET_COUNT> temp_files;
	std::vector<CompactEdge> buffer_;
	int current_bucket_;
	size_t i_, size_, node_count_;
};

}}}

#endif