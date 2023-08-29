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
#ifdef WITH_MCL
#include "../contrib/mcl/mcl.h"
#endif
#include "cascaded/cascaded.h"
#include "incremental/incremental.h"

namespace Workflow { namespace Cluster{
class ClusterRegistryStatic{
	std::map<std::string, ClusteringAlgorithm*> regMap;
public:
	ClusterRegistryStatic(){
		// To include new clustering algorithms add them into regMap
#ifdef WITH_MCL
		regMap[MCL::get_key()] = new MCL();
#endif
		regMap[::Cluster::Cascaded::get_key()] = new ::Cluster::Cascaded();		
		regMap[::Cluster::Incremental::Algo::get_key()] = new ::Cluster::Incremental::Algo();
	}
	~ClusterRegistryStatic(){
		for(auto it = regMap.begin(); it != regMap.end(); it++){
			delete it->second;
			it->second = nullptr;
		}
	}
	ClusteringAlgorithm* get(std::string key) const{
		auto ca = regMap.find(key);
		if(ca == regMap.end()){
			throw std::runtime_error("Clustering algorithm not found.");
		}
		return ca->second;
	}
	bool has(std::string key) const{
		return regMap.find(key) != regMap.end();
	}
	std::vector<std::string> getKeys() const{
		std::vector<std::string> keys;
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
	static ClusteringAlgorithm* get(std::string key){
		return reg.get(key);
	}
	static bool has(std::string key){
		return reg.has(key);
	}
	static std::vector<std::string> getKeys(){
		return reg.getKeys();
	}
};

}}
