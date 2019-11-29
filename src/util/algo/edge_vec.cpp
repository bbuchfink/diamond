#include <string>
#include "edge_vec.h"
#include "../../util/io/text_input_file.h"
#include "../../util/string/tokenizer.h"
#include "../merge_sort.h"
#include "../../basic/config.h"

using std::string;

namespace Util { namespace Algo { namespace UPGMA_MC {

EdgeVec::EdgeVec(const char *file) {
	TextInputFile in(file);
	string query, target;
	double evalue;
	while (in.getline(), !in.eof()) {
		String::Tokenizer t(in.line, "\t");
		t >> query >> target >> evalue;
		auto i = acc2idx.emplace(query, (int)acc2idx.size()), j = acc2idx.emplace(target, (int)acc2idx.size());
		if (i.second)
			idx2acc[i.first->second] = i.first->first;
		if (j.second)
			idx2acc[j.first->second] = j.first->first;
		if (i.first->second < j.first->second)
			push_back({ i.first->second, j.first->second, evalue });
	}
	merge_sort(begin(), end(), config.threads_);
}

std::string EdgeVec::print(int idx) const {
	return idx < idx2acc.size() ? idx2acc.find(idx)->second : std::to_string(idx);
}

}}}