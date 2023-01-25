#pragma once
#include "../cluster.h"

namespace Cluster { namespace Incremental {

class Algo: public ClusteringAlgorithm {
public:
	~Algo(){};
	virtual void run() override;
	std::string get_description() override;
	static std::string get_key(){
		return "incremental";
	}
};

}}
