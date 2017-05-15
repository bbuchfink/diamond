/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef COMPARE_H_
#define COMPARE_H_

#include "../basic/config.h"
#include "match_file.h"
#include "../util/compressed_stream.h"
#include "../util/seq_file_format.h"

struct Cmp_stats
{
	Cmp_stats():
		queries(0), queries1(0), queries2(0), unique1(0), unique2(0), queries1_sc(0), unique1_sc(0),
		matches1(0), matches1_hit(0), matches1_badscore(0), query_sens (0)
	{}
	size_t queries, queries1, queries2, unique1, unique2, queries1_sc, unique1_sc, matches1, matches1_hit, matches1_badscore;
	double query_sens;
};

void trim(string &s, const vector<char> &in)
{
	s.clear();
	for(size_t i=0;i<in.size();++i)
		s += in[i];
	s = s.substr(0, s.find(' '));
}

bool unique_match(const match_file::mcont::const_iterator &i, const match_file::mcont::const_iterator &begin)
{
	return i == begin || i->subject != (i-1)->subject;
}

bool consider_match(const blast_match &match)
{
	return match.n < config.max_alignments && match.bitscore >= config.min_bit_score && match.expect <= config.max_evalue;
}

void get_target_seq(vector<blast_match>::const_iterator& i,
		vector<blast_match>::const_iterator& j,
		const vector<blast_match>::const_iterator& end_i,
		const vector<blast_match>::const_iterator& end_j,
		Cmp_stats &stat)
{
	unsigned v2_matches = config.run_len != 0 ? config.run_len : 0xffffffffu;
	vector<blast_match>::const_iterator i_begin = i, j_begin = j;
	double sc_i=0,sc_j=0;
	unsigned rs_i = 0, rs_j = 0;

	while(i < end_i && i->subject == i_begin->subject) {
		sc_i = std::max(sc_i, i->bitscore);
		rs_i = std::max(rs_i, i->raw_score);
		++i;
	}

	while(j< end_j && j->subject == j_begin->subject) {
		sc_j = std::max(sc_j, j->bitscore);
		rs_j = std::max(rs_j, j->raw_score);
		++j;
	}
	if(consider_match(*i_begin)) {
		++stat.matches1;
		if(j_begin->n < v2_matches) ++stat.matches1_hit;
		double q = sc_j / sc_i;
		//if(q >= 0.95 && q <=1.05) ++matches_hit;
		//if(q < 0.95) ++matches_badscore;
		if (rs_i != rs_j) ++stat.matches1_badscore;
	}

}

void query_sens(match_file::mcont &v1,
		match_file::mcont &v2,
		Cmp_stats &stat)
{
	std::stable_sort(v1.begin(), v1.end());
	std::stable_sort(v2.begin(), v2.end());
	vector<blast_match>::const_iterator i = v1.begin(), j = v2.begin();
	size_t matches = stat.matches1, matches_hit = stat.matches1_hit;

	while(i < v1.end() && j < v2.end())
	{
		int c = i->subject.compare(j->subject);
		if(c < 0) {
			if(consider_match(*i) && unique_match(i, v1.begin())) ++stat.matches1;
			++i;
		} else if(c > 0) {
			++j;
		} else {
			get_target_seq(i, j, v1.end(), v2.end(), stat);
		}
	}
	if (stat.matches1 - matches > 0)
		stat.query_sens += double(stat.matches1_hit - matches_hit) / (stat.matches1 - matches);
}

void lone_query(match_file::mcont &v1,
		Cmp_stats &stat)
{
	vector<blast_match>::const_iterator i = v1.begin();
	while(i != v1.end()) {
		if(consider_match(*i) && unique_match(i, v1.begin()))
			++stat.matches1;
		++i;
	}
}

void print_out(Cmp_stats &stat)
{
	printf("queries=%zu queries(1)=%zu queries(2)=%zu\n", stat.queries, stat.queries1, stat.queries2);
	printf("unique(1)=%zu unique(2)=%zu\n", stat.unique1, stat.unique2);
	printf("queries(1)>sc=%zu unique(1)>sc=%zu hit(2)=%zu (%1.f%%)\n", stat.queries1_sc, stat.unique1_sc, stat.queries1_sc-stat.unique1_sc, (double)(stat.queries1_sc-stat.unique1_sc)*100/stat.queries1_sc);
	printf("matches(1)>sc=%zu hit(2)=%zu (%.1lf%%) bad score=%zu\n", stat.matches1, stat.matches1_hit, double(stat.matches1_hit)*100/stat.matches1, stat.matches1_badscore);
	printf("query_sens=%.1lf\n", stat.query_sens * 100 / stat.queries1_sc);
	printf("\n");
}

void compare()
{
	typedef blast_tab_format_with_rawscore Format1;
	typedef blast_tab_format Format2;

	vector<char> id;
	vector<Letter> seq;

	Input_stream seqStream(config.query_file);
	match_file file1 (config.match_file1.c_str());
	match_file::mcont v1;
	file1.get_read(v1, Format1());
	match_file file2 (config.match_file2.c_str());
	match_file::mcont v2;
	file2.get_read(v2, Format2());
	bool do_out = config.output_file.length() > 0;
	FILE *out = 0;
	if(do_out) out = fopen(config.output_file.c_str(), "wt");

	Cmp_stats stat;
	string q;
	size_t read = 0;
	FASTA_format format;
	while(format.get_seq(id, seq, seqStream)) {
		trim(q, id);
		++stat.queries;
		//printf("%lu ", queries);
		if (stat.queries % 1000 == 0) {
			printf("n = %zu\n", stat.queries);
			print_out(stat);
		}

		bool have1 = false, have2 = false, have1_sc = false;
		if(v1.size() > 0 && q == v1[0].query) {
			++stat.queries1;
			if(v1[0].bitscore >= config.min_bit_score && v1[0].expect <= config.max_evalue) {
				have1_sc = true;
				++stat.queries1_sc;
			}
			have1 = true;
		}

		if(v2.size() > 0 && q == v2[0].query) {
			++stat.queries2;
			have2 = true;
		} else if(have1_sc) {
			++stat.unique1_sc;
		}

		if(have1_sc && !have2)
			lone_query(v1, stat);
		else if(have1_sc && have2)
			query_sens(v1, v2, stat);

		if(have1 && !have2) {
			++stat.unique1;
			if(do_out) fprintf(out, "1 %s\n", q.c_str());
		} else if(have2 && !have1) {
			++stat.unique2;
			if(do_out) fprintf(out, "2 %s\n", q.c_str());
		}
		if(have1) file1.get_read(v1, Format1());
		if(have2) file2.get_read(v2, Format2());

		++read;
	}

	print_out(stat);

	if(do_out) fclose(out);
}

#endif /* COMPARE_H_ */
