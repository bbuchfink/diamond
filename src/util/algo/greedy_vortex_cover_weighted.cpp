#include <algorithm>
#include "algo.h"

using namespace std;

namespace Util { namespace Algo {

int fix_centroid(vector<int> &centroid, int i) {
	if (centroid[i] == i)
		return i;
	else
		return centroid[i] = fix_centroid(centroid, centroid[i]);
}

vector<int> greedy_vortex_cover_weighted(vector<Edge> &edges, int vortex_count) {
	sort(edges.begin(), edges.end());
	vector<int> centroid, neighbors(vortex_count);
	centroid.reserve(vortex_count);
	for (int i = 0; i < vortex_count; ++i)
		centroid.push_back(i);
	for (Edge e : edges) {
		++neighbors[e.v1];
		++neighbors[e.v2];
	}
	for (Edge e : edges) {
		if (centroid[e.v1] == e.v1 && centroid[e.v2] == e.v2) {
			if (neighbors[e.v1] >= neighbors[e.v2])
				centroid[e.v2] = e.v1;
			else
				centroid[e.v1] = e.v2;
		}
	}
	for (int i = 0; i < vortex_count; ++i)
		fix_centroid(centroid, i);
	return centroid;
}

}}