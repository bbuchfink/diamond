#pragma once
#include <string>
class ClusteringAlgorithm {
public:
	virtual void run() = 0;
	static std::string get_key();
};
