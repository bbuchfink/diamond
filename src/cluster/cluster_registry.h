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

namespace Workflow { namespace Cluster{
class ClusterRegistryStatic{
	std::map<std::string, ClusteringAlgorithm*> regMap;
public:
	ClusterRegistryStatic(){
		// To include new clustering algorithms add them into regMap
		regMap[MultiStep::get_key()] = new MultiStep();
		regMap[MCL::get_key()] = new MCL();
	}
	~ClusterRegistryStatic(){
		for(auto it = regMap.begin(); it != regMap.end(); it++){
			delete it->second;
			it->second = nullptr;
		}
	}
	ClusteringAlgorithm* get(string key) const{
		auto ca = regMap.find(key);
		if(ca == regMap.end()){
			throw std::runtime_error("Clustering algorithm not found.");
		}
		return ca->second;
	}
	bool has(string key) const{
		return regMap.find(key) != regMap.end();
	}
	vector<string> getKeys() const{
		vector<string> keys;
		for(auto it = regMap.begin(); it != regMap.end(); it++){
			keys.push_back(it->first);
		}
		return keys;
	}

};

class ClusterRegistry{
private:
	static const ClusterRegistryStatic reg;
public:
	static ClusteringAlgorithm* get(string key){
		return reg.get(key);
	}
	static bool has(string key){
		return reg.has(key);
	}
	static vector<string> getKeys(){
		return reg.getKeys();
	}
};

}}
