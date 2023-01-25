#include "cascaded.h"

using std::string;
using std::vector;

namespace Cluster {

vector<string> cluster_steps(double approx_id) {
	if (!config.cluster_steps.empty())
		return config.cluster_steps;
	vector<string> v = { "faster_lin", "fast" };
	if (approx_id < 90)
		v.push_back("default");
	if (approx_id < 50)
		v.push_back("more-sensitive");
	return v;
}

}