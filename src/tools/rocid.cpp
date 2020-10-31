#include <array>
#include <unordered_map>
#include "../basic/config.h"
#include "../util/io/text_input_file.h"
#include "../util/string/tokenizer.h"
#include "../util/log_stream.h"

using std::cout;
using std::endl;
using std::array;
using std::unordered_map;

static string query_aln;
static array<int, 10> totals, counts;
unordered_map<string, uint8_t> acc2id;
static size_t unmapped_query = 0;

static void fetch_map(TextInputFile& map_in, const string& query) {
	string q, target;
	float id;
	acc2id.clear();
	std::fill(totals.begin(), totals.end(), 0);
	while (map_in.getline(), !map_in.eof()) {
		Util::String::Tokenizer(map_in.line, "\t") >> q >> target >> id;
		if (q != query) {
			if (q < query)
				continue;
			else {
				map_in.putback_line();
				return;
			}
		}
		int bin = std::min((int)(id*100.0) / 10, 9);
		acc2id[target] = (uint8_t)bin;
		++totals[bin];
	}
}

static void print() {
	if (unmapped_query)
		return;
	cout << query_aln;
	for (int i = 0; i < 10; ++i)
		cout << '\t' << (totals[i] ? ((double)counts[i] / (double)totals[i]) : -1.0);
	cout << endl;
}

void roc_id() {
	TextInputFile in(config.query_file);
	string query, target;
	size_t n = 0, queries = 0, unmapped = 0, hits = 0;
	std::fill(counts.begin(), counts.end(), 0);

	TextInputFile map_in(config.family_map);
	
	while (in.getline(), !in.eof()) {
		if (in.line.empty())
			break;
		Util::String::Tokenizer(in.line, "\t") >> query >> target;
		++hits;
		if (query != query_aln) {
			print();
			fetch_map(map_in, query);
			query_aln = query;
			++queries;
			std::fill(counts.begin(), counts.end(), 0);
			unmapped_query = 0;
			if (queries % 1000 == 0)
				message_stream << queries << ' ' << hits << ' ' << unmapped << endl;
		}

		auto it = acc2id.find(target);
		if (it == acc2id.end()) {
			++unmapped_query;
			++unmapped;
			continue;
		}
		++counts[(int)it->second];
	}
	print();

	in.close();
	map_in.close();
	message_stream << "Queries = " << queries << endl;

}