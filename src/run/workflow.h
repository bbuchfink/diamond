#ifndef RUN_WORKFLOW_H_
#define RUN_WORKFLOW_H_

#include <vector>

struct DatabaseFile;
struct Consumer;
struct TextInputFile;

namespace Workflow { 
namespace Search {

struct Options {
	Options():
		self(false),
		db(nullptr),
		consumer(nullptr),
		query_file(nullptr),
		db_filter(nullptr)
	{}
	bool self;
	DatabaseFile *db;
	Consumer *consumer;
	TextInputFile *query_file;
	const std::vector<bool> *db_filter;
};

void run(const Options &options);

}
}

#endif
