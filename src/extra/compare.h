#ifndef COMPARE_H_
#define COMPARE_H_

#include "../basic/config.h"
#include "match_file.h"
#include "../util/compressed_stream.h"
#include "../util/seq_file_format.h"

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
		size_t &matches,
		size_t &matches_hit,
		size_t &matches_badscore)
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
		++matches;
		if(j_begin->n < v2_matches) ++matches_hit;
		double q = sc_j / sc_i;
		//if(q >= 0.95 && q <=1.05) ++matches_hit;
		//if(q < 0.95) ++matches_badscore;
		if (rs_i != rs_j) ++matches_badscore;
	}

}

void query_sens(match_file::mcont &v1,
		match_file::mcont &v2,
		size_t &matches,
		size_t &matches_hit,
		size_t &matches_badscore)
{
	std::stable_sort(v1.begin(), v1.end());
	std::stable_sort(v2.begin(), v2.end());
	vector<blast_match>::const_iterator i = v1.begin(), j = v2.begin();

	while(i < v1.end() && j < v2.end())
	{
		int c = i->subject.compare(j->subject);
		if(c < 0) {
			if(consider_match(*i) && unique_match(i, v1.begin())) ++matches;
			++i;
		} else if(c > 0) {
			++j;
		} else {
			get_target_seq(i, j, v1.end(), v2.end(), matches, matches_hit, matches_badscore);
		}
	}
}

void lone_query(match_file::mcont &v1,
		size_t &matches)
{
	vector<blast_match>::const_iterator i = v1.begin();
	while(i != v1.end()) {
		if(consider_match(*i) && unique_match(i, v1.begin()))
			++matches;
		++i;
	}
}

void print_out(size_t queries,size_t queries1, size_t queries2, size_t unique1, size_t unique2, size_t queries1_sc, size_t unique1_sc, size_t matches1, size_t matches1_hit, size_t matches1_badscore)
{
	printf("queries=%zu queries(1)=%zu queries(2)=%zu\n", queries, queries1, queries2);
	printf("unique(1)=%zu unique(2)=%zu\n", unique1, unique2);
	printf("queries(1)>sc=%zu unique(1)>sc=%zu hit(2)=%zu (%1.f%%)\n", queries1_sc, unique1_sc, queries1_sc-unique1_sc, (double)(queries1_sc-unique1_sc)*100/queries1_sc);
	printf("matches(1)>sc=%zu hit(2)=%zu (%.1lf%%) bad score=%zu\n", matches1, matches1_hit, double(matches1_hit)*100/matches1, matches1_badscore);
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

	size_t queries (0), queries1 (0), queries2 (0), unique1 (0), unique2 (0), queries1_sc (0), unique1_sc (0),
			matches1 (0), matches1_hit (0), matches1_badscore (0);
	string q;
	size_t read = 0;
	FASTA_format format;
	while(format.get_seq(id, seq, seqStream)) {
		trim(q, id);
		++queries;
		//printf("%lu ", queries);
		if(queries % 1000 == 0) {
			printf("n = %zu\n", queries);
			print_out(queries,queries1, queries2, unique1, unique2, queries1_sc, unique1_sc, matches1, matches1_hit, matches1_badscore);
		}

		bool have1 = false, have2 = false, have1_sc = false;
		if(v1.size() > 0 && q == v1[0].query) {
			++queries1;
			if(v1[0].bitscore >= config.min_bit_score && v1[0].expect <= config.max_evalue) {
				have1_sc = true;
				++queries1_sc;
			}
			have1 = true;
		}

		if(v2.size() > 0 && q == v2[0].query) {
			++queries2;
			have2 = true;
		} else if(have1_sc) {
			++unique1_sc;
		}

		if(have1_sc && !have2)
			lone_query(v1, matches1);
		else if(have1_sc && have2)
			query_sens(v1, v2, matches1, matches1_hit, matches1_badscore);

		if(have1 && !have2) {
			++unique1;
			if(do_out) fprintf(out, "1 %s\n", q.c_str());
		} else if(have2 && !have1) {
			++unique2;
			if(do_out) fprintf(out, "2 %s\n", q.c_str());
		}
		if(have1) file1.get_read(v1, Format1());
		if(have2) file2.get_read(v2, Format2());

		++read;
	}

	print_out(queries,queries1, queries2, unique1, unique2, queries1_sc, unique1_sc, matches1, matches1_hit, matches1_badscore);

	if(do_out) fclose(out);
}

#endif /* COMPARE_H_ */
