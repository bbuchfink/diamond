#include <iostream>
#include <array>
#include <set>
#include "../util/io/text_input_file.h"
#include "../basic/config.h"
#include "../util/intrin.h"

using std::set;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::array;

static const int W = 7, L = 16, N = 64;

static array<size_t, 1 << L> counts;
static set<int> exclude;

static int pattern(int p) {
	while ((p & 1) == 0)
		p >>= 1;
	return p;
}

static void process_pattern(const string& s, set<int>& patterns) {
	int p = 0;
	for (size_t i = 0; i < s.length(); ++i) {
		p <<= 1;
		if (s[i] == '1')
			p |= 1;
		p &= (1 << L) - 1;
		if (popcount32((unsigned)p) == W) {
			patterns.insert(pattern(p));
		}
	}
}

static void process_aln(vector<string>::const_iterator begin, vector<string>::const_iterator end) {
	set<int> patterns;
	for (auto i = begin; i < end; ++i)
		process_pattern(*i, patterns);

	vector<int> v;
	std::set_intersection(exclude.begin(), exclude.end(), patterns.begin(), patterns.end(), std::back_inserter(v));
	if(v.empty())
		for (int p : patterns)
			++counts[p];
}

void find_shapes() {
	int hit_total = 0;
	for (int i = 0; i < N; ++i) {
		counts.fill(0);
		TextInputFile in(config.single_query_file());
		size_t n = 0;
		while (in.getline(), !in.eof() || !in.line.empty()) {
			auto v = tokenize(in.line.c_str(), "\t");
			process_aln(v.begin(), v.end());
			++n;
		}
		in.close();
		int p_max = 0;
		int c_max = 0;
		for (int p = 0; p < (int)counts.size(); ++p)
			if (counts[p] > p_max) {
				p_max = p;
				c_max = counts[p];
			}
		cout << "Alignments: " << n << endl;
		cout << "Pattern: " << p_max << endl;
		cout << "Hit: " << c_max << endl;
		hit_total += c_max;
		cout << "Total hit: " << hit_total << endl;
		exclude.insert(p_max);
	}
}