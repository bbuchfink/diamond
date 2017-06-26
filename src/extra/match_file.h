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

#ifndef MATCH_FILE_H_
#define MATCH_FILE_H_

#include <stdio.h>
#include <vector>
#include "../util/binary_file.h"
#include "blast_record.h"
#include "../util/util.h"
#include "../basic/reduction.h"
#include "../basic/shape_config.h"

struct file_parse_exception : public std::runtime_error
{
	file_parse_exception(size_t l):
		std::runtime_error (string("File parse error line ") + to_string((unsigned)l) + "\n")
	{ }
};

struct blast_format { };
struct blast_tab_format { };
struct blast_tab_format_with_rawscore {};

using std::vector;

inline void alline(const char *buffer, char *dest)
{
	buffer += 7;
	const char *begin = strstr(buffer, " ") + 1;
	while(*begin == ' ')
		++begin;
	const char *end = strstr(begin, " ");
	strncpy(dest, begin, end - begin);
	dest[end-begin] = 0;
}

inline bool seq_letter(char x)
{
	return x >= 'A' && x <= 'Z';
}

/*int64_t match_scores[32];
int64_t match_counts[32];*/

inline unsigned char_match(char query, char subject, const Reduction &red)
{
	Letter ql = query;
	Letter sl = subject;
	if(red(ql) == red(sl)) {
		//int sc = score_matrix::get().letter_score(ql, sl);
		//match_scores[(int)ql] += sc;
		//++match_counts[(int)ql];
		//match_scores[(int)sl] += sc;
		//++match_counts[(int)sl];
		return 1;
	} else
		return 0;
}

inline void get_match(unsigned &mask, const char *queryl, const char *subjectl, bool &hit, unsigned &len, unsigned &rid, unsigned &current_len, unsigned &ungapped_len, Letter* q, Letter* s, bool &stop)
{
	unsigned l ((unsigned)strlen(queryl));
	for(unsigned i=0;i<l;++i) {
		mask <<= 1;
		if(seq_letter(queryl[i]) && seq_letter(subjectl[i])) {
			unsigned x = char_match(queryl[i], subjectl[i], Reduction::reduction);
			rid += x;
			++len;
			++current_len;
			ungapped_len = std::max(ungapped_len, current_len);
			mask |= x;
		} else if(queryl[i] == '-' || subjectl[i] == '-') {
			mask = 0;
			current_len = 0;
		} else if(queryl[i] == '*') {
			stop = true;
		}
		*(q++) = queryl[i];
		*(s++) = subjectl[i];
		for(unsigned j=0;j<shapes.count();++j) {
			if(((mask & shapes[j].rev_mask_) == shapes[j].rev_mask_))
				hit = true;
		}
	}
}

class match_file : public Input_stream
{

public:

	typedef std::vector<blast_match> mcont;

	match_file(const char* fileName) : Input_stream(fileName), currentQueryCount(0), queryCount(0), matchCount(0)
	{
		currentQuery[0] = 0;
		currentSubject[0] = 0;
		memset(subst_p, 0, sizeof(subst_p));
		memset(subst_n, 0, sizeof(subst_n));
	}

	void set_subst(const char* q, const char *s)
	{
		while (*q) {
			if (*q != '-' && *s != '-' && *q != *s) {
				const Letter lq = value_traits.from_char(*q), ls = value_traits.from_char(*s);
				if (lq < 20 && ls < 20 && lq != ls) {
					++subst_p[(size_t)lq][(size_t)ls];
					++subst_n[(size_t)lq];
				}
			}
			++q;
			++s;
		}
	}

	void get_subst() const
	{
		for (unsigned i = 0; i < 20; ++i) {
			cout << "{";
			for (unsigned j = 0; j < 20; ++j)
				cout << (double)subst_p[i][j] / subst_n[i] << ',';
			cout << "}," << endl;
		}
	}

	bool get(blast_match &record, const blast_format&)
	{

		enum State { begin = 0, end = 1, queryStart = 2, subjectStart = 3, matchStart = 4, queryLine = 5, subjectLine = 6, separator = 7, between = 8, haveid = 9 } state = begin;
		unsigned expect_i, id1, id2;
		char queryl[128], subjectl[128];
		Letter queryseq[4096], subjectseq[4096];
		Letter *q = queryseq, *s = subjectseq;
		unsigned match_mask = 0, current_len = 0;
		record.hit = false;
		record.ungapped_len = 0;
		record.rid = 0;
		record.len = 0;
		record.stop = false;

		while(state != end && (this->getline(), !this->eof())) {

			//printf("%lu %u %s\n", line_count, state, line.c_str());

			if(!strncmp(line.c_str(), "Query= ", 7)) {
				if(state == begin || state == queryStart) {
					state = queryStart;
					if(sscanf(line.c_str(), "Query= %s", currentQuery) != 1)
						throw file_parse_exception(this->line_count);
				} else
					throw file_parse_exception(this->line_count);
			} else if(line.c_str()[0] == '>') {
				if(state == begin || state == queryStart) {
					state = subjectStart;
					if(sscanf(line.c_str(), "> %s", currentSubject) != 1)
						throw file_parse_exception(this->line_count);
				}
				else if (state == separator) {
					this->putback_line();
					state = end;
				}
				else
					throw file_parse_exception(this->line_count);
			} else if(sscanf(line.c_str(), " Score = %lf bits (%u),  Expect = %lf", &record.bitscore, &record.raw_score, &record.expect) == 3
					|| sscanf(line.c_str(), " Score = %lf bits (%u),  Expect(%i) = %lf", &record.bitscore, &record.raw_score, &expect_i, &record.expect) == 4) {
				if(state == subjectStart || state == begin) {
					record.query = currentQuery;
					record.subject = currentSubject;
					state = matchStart;
				}
				else if (state == separator) {
					this->putback_line();
					state = end;
				}
				else
					throw file_parse_exception(this->line_count);
			} else if(sscanf(line.c_str(), " Identities = %u/%u (%lf%%)", &id1, &id2, &record.id) == 3) {
				if(state == matchStart) {
					state = haveid;
				} else
					throw file_parse_exception(this->line_count);
			} else if(!strncmp(line.c_str(), "Query", 5)) {
				if(state == haveid || state == separator) {
					state = queryLine;
					alline(line.c_str(), queryl);
				} else
					throw file_parse_exception(this->line_count);
			} else if(!strncmp(line.c_str(), "Sbjct", 5)) {
				if(state == between) {
					state = subjectLine;
					alline(line.c_str(), subjectl);
					set_subst(queryl, subjectl);
					/*if(strlen(subjectl) != strlen(queryl))
						throw file_parse_exception(lineNumber);*/
					//get_match(match_mask, queryl, subjectl, record.hit, record.len, record.rid, current_len, record.ungapped_len, q, s, record.stop);
				} else
					throw file_parse_exception(this->line_count);
			} else if(state == queryLine && line.c_str()[0] == ' ') {
				state = between;
			} else if((state == subjectLine || state == separator) && line.empty()) {
				//printf("%i\n",(int)buffer[0]);
				if(state == subjectLine)
					state = separator;
				else if(state == separator)
					state = end;
				else
					throw file_parse_exception(this->line_count);
			} else {
				if(state != begin && state != queryStart && state != haveid && state != subjectStart)
					throw file_parse_exception(this->line_count);
			}

		}

		if(this->eof())
			return false;
		if(state != end)
			throw file_parse_exception(this->line_count);
		++currentQueryCount;
		++matchCount;
		return true;

	}

	bool get(blast_match &match, const blast_tab_format&)
	{
		char query[nameBufferSize], s2[nameBufferSize];
		double f1;
		unsigned long long i1,i2,i3,i4,i5,i6,i7;
		this->getline();
		if(!this->eof()) {
			while (this->line[0] == '#')
				this->getline();
			if(sscanf(line.c_str(), "%s%s%lf%llu%llu%llu%llu%llu%llu%llu%lf%lf", query, s2, &f1, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &match.expect, &match.bitscore) == 12) {
				//|| sscanf(line.c_str(), "%s%s%s%lf%llu%llu%llu%llu%llu%llu%llu%lf%lf", query, query2, s2, &f1, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &match.expect, &match.bitscore) == 13) {
				/*++matchCount;
				if(!strcmp(query, currentQuery)) {
					++currentQueryCount;
				} else {
					++queryCount;
					currentQueryCount = 1;
					strcpy(currentQuery, query);
				}*/
				match.query = query;
				match.subject = s2;
				return true;
			} else {
				printf("%s\n", line.c_str());
				throw file_parse_exception(this->line_count);
			}
		} else {
			match.set_empty();
			return false;
		}
	}

	bool get(blast_match &match, const blast_tab_format_with_rawscore&)
	{
		char query[nameBufferSize], s2[nameBufferSize];
		double f1;
		unsigned long long i1, i2, i3, i4, i5, i6, i7;
		this->getline();
		if (!this->eof()) {
			while (this->line[0] == '#')
				this->getline();
			if (sscanf(line.c_str(), "%s%s%lf%llu%llu%llu%llu%llu%llu%llu%lf%lf%u", query, s2, &f1, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &match.expect, &match.bitscore, &match.raw_score) == 13) {
				//|| sscanf(line.c_str(), "%s%s%s%lf%llu%llu%llu%llu%llu%llu%llu%lf%lf", query, query2, s2, &f1, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &match.expect, &match.bitscore) == 13) {
				/*++matchCount;
				if(!strcmp(query, currentQuery)) {
				++currentQueryCount;
				} else {
				++queryCount;
				currentQueryCount = 1;
				strcpy(currentQuery, query);
				}*/
				match.query = query;
				match.subject = s2;
				return true;
			}
			else {
				printf("%s\n", line.c_str());
				throw file_parse_exception(this->line_count);
			}
		}
		else {
			match.set_empty();
			return false;
		}
	}

	template<typename _format>
	bool get_read(mcont &v, const _format& format)
	{
		unsigned n = 0;
		v.clear();
		blast_match match;
		if(!save.empty()) {
			save.n = n++;
			v.push_back(save);
			save.set_empty();
		} else {
			if(!get(match, format))
				return false;
			match.n = n++;
			v.push_back(match);
		}

		while(get(match, format) && match.query == v[0].query) {
			if(match.subject != v.back().subject)
				n++;
			match.n = n;
			v.push_back(match);
		}
		if(!match.empty() && match.query != v[0].query)
			save = match;
		return !v.empty();
	}

	unsigned getTotalQueries() const
	{
		return queryCount;
	}

	unsigned getCurrentQueryCount() const
	{
		return currentQueryCount;
	}

protected:
	static const size_t nameBufferSize = 4096;
	unsigned currentQueryCount, queryCount, matchCount;
	char currentQuery[nameBufferSize], currentSubject[nameBufferSize];
	blast_match save;
	size_t subst_p[20][20];
	size_t subst_n[20];

};

#endif /* MATCH_FILE_H_ */
