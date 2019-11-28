#ifndef EDGE_VEC_H_
#define EDGE_VEC_H_

#include <unordered_map>
#include <vector>
#include <string>

namespace Util { namespace Algo { namespace UPGMA_MC {

struct CompactEdge {
	bool operator<(const CompactEdge &e) const {
		return d < e.d;
	}
	int n1, n2;
	float d;
};

struct EdgeVec : public std::vector<CompactEdge> {
	EdgeVec(const char *file);
	std::unordered_map<std::string, int> acc2idx;
};

}}}

#endif