#include <queue>
#include "file.h"

using std::priority_queue;
using std::pair;
using std::vector;
using std::greater;

namespace Util { namespace Tsv {

File* merge(std::vector<File*>::iterator begin, std::vector<File*>::iterator end, int column) {
	using pair = pair<int64_t, int64_t>;
	File* out = new File((*begin)->schema(), "", Flags::TEMP);
	priority_queue<pair, vector<pair>, greater<pair>> queue;
	vector<Table> tables;
	tables.reserve(end - begin);
	for (auto it = begin; it != end; ++it) {
		tables.push_back((*it)->read_record());
		if (tables.back().empty())
			continue;
		queue.emplace(tables.back().front().get<int64_t>(column), it - begin);
	}
	while (!queue.empty()) {
		const int64_t f = queue.top().second;
		out->write(tables[f].front());
		queue.pop();
		tables[f] = begin[f]->read_record();
		if (!tables[f].empty())
			queue.emplace(tables[f].front().get<int64_t>(column), f);
	}
	return out;
}

}}