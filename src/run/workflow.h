#ifndef RUN_WORKFLOW_H_
#define RUN_WORKFLOW_H_

#include <vector>

struct DatabaseFile;
struct Consumer;

namespace Workflow { 

namespace Search {

struct Options {
	Options():
		self(false),
		db(nullptr),
		consumer(nullptr),
		db_filter(nullptr)
	{}
	bool self;
	DatabaseFile *db;
	Consumer *consumer;
	const std::vector<bool> *db_filter;
};

void run(const Options &options);

}

namespace Cluster {

void run();

}

}

#endif