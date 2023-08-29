#include <iostream>
#include <array>
#include <set>
#include "../util/io/text_input_file.h"
#include "../basic/config.h"
#include "../util/intrin.h"
#include "../util/util.h"

using std::set;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::array;

using Pattern = uint32_t;

static const int W = 7, L = 16, N = 64, T = 3;

static array<size_t, 1 << L> counts;
static set<Pattern> exclude;

static void process_window(int p, set<Pattern>& patterns) {
	int n = 0;
	array<int, L> idx;
	for (int i = 1; i < L; ++i)
		if (p & (1 << i))
			idx[n++] = i;

	string bitmask(W - 1, 1);
	bitmask.resize(n, 0);

	do {
		int q = 1;
		for (int i = 0; i < n; ++i)
			if (bitmask[i]) q |= 1 << idx[i];
		patterns.insert(q);
	} while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}

static void process_pattern(const string& s, set<Pattern>& patterns) {
	int p = 0;
	for (size_t i = 0; i < s.length(); ++i) {
		p <<= 1;
		p &= (1 << L) - 1;
		if (s[i] == '1') {
			p |= 1;			
			const int n = popcount32((unsigned)p);
			if (n >= W && n <= W + T) {
				process_window(p, patterns);
			}
		}
	}
}

static bool is_excluded(vector<string>::const_iterator begin, vector<string>::const_iterator end) {
	if (exclude.empty())
		return false;
	for (auto it = begin; it < end; ++it) {
		int p = 0;
		for (size_t i = 0; i < it->length(); ++i) {
			p <<= 1;
			if ((*it)[i] == '1')
				p |= 1;
			p &= (1 << L) - 1;
			const int n = popcount32((unsigned)p);
			if (n >= W)
				for (int e : exclude)
					if ((p & e) == e)
						return true;
		}
	}
	return false;
}

static bool is_id(const string& s) {
	size_t n = 0;
	for (char c : s)
		if (c == '1')
			++n;
	return n == s.length();
}

static void process_aln(vector<string>::const_iterator begin, vector<string>::const_iterator end) {
	set<Pattern> patterns;
	for (auto i = begin; i < end; ++i)
		process_pattern(*i, patterns);

	vector<int> v;
	std::set_intersection(exclude.begin(), exclude.end(), patterns.begin(), patterns.end(), std::back_inserter(v));
	if(v.empty())
		for (int p : patterns)
			++counts[p];
}

static string as_string(Pattern p) {
	string s;
	for (int i = 0; i < L; ++i)
		if (p & (Pattern(1) << i))
			s = '1' + s;
		else
			s = '0' + s;
	while (s.front() == '0')
		s.erase(s.begin());
	return s;
}

static void print_all() {
	for (int p : exclude) {
		cout << as_string(p) << endl;
	}
}

void find_shapes() {
	for (int i = 0; i < N; ++i) {
		counts.fill(0);
		TextInputFile in(config.single_query_file());
		size_t n = 0, nt = 0;
		while (in.getline(), !in.eof() || !in.line.empty()) {
			auto v = tokenize(in.line.c_str(), "\t");
			if ((v.size() > 1 || !is_id(v.front()))) {
				if (!is_excluded(v.begin(), v.end())) {
					process_aln(v.begin(), v.end());
					++n;
				}
				++nt;
			}
			if ((nt % 100000) == 0)
				cout << "Processed = " << nt << endl;
		}
		in.close();
		Pattern p_max = 0;
		size_t c_max = 0;
		for (Pattern p = 0; p < (Pattern)counts.size(); ++p)
			if (counts[p] > c_max) {
				p_max = p;
				c_max = counts[p];
			}
		cout << "Alignments: " << n << " / " << nt << endl;
		cout << "Pattern: " << as_string(p_max) << endl;
		cout << "Hit: " << c_max << endl;
		exclude.insert(p_max);
		print_all();
	}	
}