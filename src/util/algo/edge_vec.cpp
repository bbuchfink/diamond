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
		auto i = acc2idx.emplace(query, (int)acc2idx.size()).first, j = acc2idx.emplace(target, (int)acc2idx.size()).first;
		if (i->second < j->second)
			push_back({ i->second, j->second, evalue });
	}
	merge_sort(begin(), end(), config.threads_);
}

}}}