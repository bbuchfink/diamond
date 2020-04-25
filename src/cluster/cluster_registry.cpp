#include "cluster_registry.h"

using namespace std;

namespace Workflow { namespace Cluster{
	map<string, ClusteringAlgorithm*> ClusterRegistry::regMap;
	MCL ClusterRegistry::mcl;
	MultiStep ClusterRegistry::multiStep;
	ClusterRegistry::StaticConstructor ClusterRegistry::_staticConstructor;
}}
