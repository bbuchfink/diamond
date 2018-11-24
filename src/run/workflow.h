#ifndef RUN_WORKFLOW_H_
#define RUN_WORKFLOW_H_

struct DatabaseFile;

namespace Workflow { 

namespace Search {

struct Options {
	Options():
		self(false),
		db(nullptr)
	{}
	bool self;
	DatabaseFile *db;
};

void run(const Options &options);

}

}

#endif