#pragma once
#include "cluster.h"
#include "mcl.h"
#include "multi_step_cluster.h"

using namespace std;

namespace Workflow { namespace Cluster{
class ClusterRegistry{
private:
	ClusterRegistry(){};
	static MCL mcl;
	static MultiStep multiStep;
public:
	static map<string, ClusteringAlgorithm*> regMap;
	static ClusteringAlgorithm* get(string key){
		auto ca = ClusterRegistry::regMap.find(key);
		if(ca == ClusterRegistry::regMap.end()){
			throw std::runtime_error("Clustering algorithm not found.");
		}
		return ca->second;
	}
	static bool has(string key){
		return ClusterRegistry::regMap.find(key) != ClusterRegistry::regMap.end();
	}
	static vector<string> getKeys(){
		auto it = regMap.begin();
		vector<string> keys;
		while(it != regMap.end()){
			keys.push_back(it->first);
			it++;
		}
		return keys;
	}
	static struct StaticConstructor {
		StaticConstructor() {
			regMap.emplace(multiStep.get_key(), &multiStep);
			regMap.emplace(mcl.get_key(), &mcl);
		}
	} _staticConstructor;
};

}}
