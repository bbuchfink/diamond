/****
DIAMOND protein aligner
Copyright (C) 2020 QIAGEN A/S (Aarhus, Denmark)
Code developed by Patrick Ettenhuber <patrick.ettenhuber@qiagen.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#pragma once
#include "cluster.h"
#include "mcl.h"
#include "multi_step_cluster.h"

using namespace std;

namespace Workflow { namespace Cluster{
class ClusterRegistry{
private:
	ClusterRegistry(){};
	// To include new clustering algorithms add the instantiation here and in the cluster_registry.cpp file. Then add it to the StaticConstructor below
	static MCL mcl;
	static MultiStep multiStep;
public:
	static map<string, ClusteringAlgorithm*> regMap;
	static ClusteringAlgorithm* get(string key){
		map<string, ClusteringAlgorithm*>::iterator ca = ClusterRegistry::regMap.find(key);
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
			// Add any new clustering algorithm here
			regMap.emplace(multiStep.get_key(), &multiStep);
			regMap.emplace(mcl.get_key(), &mcl);
		}
	} _staticConstructor;
};

}}
