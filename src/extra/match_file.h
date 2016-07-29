#ifndef MATCH_FILE_H_
#define MATCH_FILE_H_

#include <vector>
#include "text_file.h"
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

using std::vector;

void alline(const char *buffer, char *dest)
{
	buffer += 7;
	const char *begin = strstr(buffer, " ") + 1;
	while(*begin == ' ')
		++begin;
	const char *end = strstr(begin, " ");
	strncpy(dest, begin, end - begin);
	dest[end-begin] = 0;
}

bool seq_letter(char x)
{
	return x >= 'A' && x <= 'Z';
}

int64_t match_scores[32];
int64_t match_counts[32];

unsigned char_match(char query, char subject, const Reduction &red)
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

void get_match(unsigned &mask, const char *queryl, const char *subjectl, bool &hit, unsigned &len, unsigned &rid, unsigned &current_len, unsigned &ungapped_len, Letter* q, Letter* s, bool &stop)
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
			if(match_shape_mask(mask, shapes[j].rev_mask_) && shapes[j].is_low_freq_rev(q))
				hit = true;
		}
	}
}

class match_file : public text_file
{

public:

	typedef std::vector<blast_match> mcont;

	match_file(const char* fileName) : text_file(fileName), currentQueryCount(0), queryCount(0), matchCount(0)
	{
		currentQuery[0] = 0;
		currentSubject[0] = 0;
	}

	bool get(blast_match &record, const blast_format&)
	{

		enum State { begin=0, end=1, queryStart=2, subjectStart=3, matchStart=4, queryLine=5, subjectLine=6, separator=7, between=8, haveid=9 } state = begin;
		unsigned rawScore, expect_i, id1, id2;
		char queryl[128], subjectl[128];
		Letter queryseq[4096], subjectseq[4096];
		Letter *q = queryseq, *s = subjectseq;
		unsigned match_mask = 0, current_len = 0;
		record.hit = false;
		record.ungapped_len = 0;
		record.rid = 0;
		record.len = 0;
		record.stop = false;

		while(state != end && readLine(buffer)) {

			//printf("%lu %u %s", lineNumber, state, buffer);

			if(!strncmp(buffer, "Query= ", 7)) {
				if(state == begin || state == queryStart) {
					state = queryStart;
					if(sscanf(buffer, "Query= %s", currentQuery) != 1)
						throw file_parse_exception(lineNumber);
				} else
					throw file_parse_exception(lineNumber);
			} else if(buffer[0] == '>') {
				if(state == begin || state == queryStart) {
					state = subjectStart;
					if(sscanf(buffer, "> %s", currentSubject) != 1)
						throw file_parse_exception(lineNumber);
				} else
					throw file_parse_exception(lineNumber);
			} else if(sscanf(buffer, " Score = %lf bits (%u),  Expect = %lf", &record.bitscore, &rawScore, &record.expect) == 3
					|| sscanf(buffer, " Score = %lf bits (%u),  Expect(%i) = %lf", &record.bitscore, &rawScore, &expect_i, &record.expect) == 4) {
				if(state == subjectStart || state == begin) {
					record.query = currentQuery;
					record.subject = currentSubject;
					state = matchStart;
				} else
					throw file_parse_exception(lineNumber);
			} else if(sscanf(buffer, " Identities = %u/%u (%lf%%)", &id1, &id2, &record.id) == 3) {
				if(state == matchStart) {
					state = haveid;
				} else
					throw file_parse_exception(lineNumber);
			} else if(!strncmp(buffer, "Query", 5)) {
				if(state == haveid || state == separator) {
					state = queryLine;
					alline(buffer, queryl);
				} else
					throw file_parse_exception(lineNumber);
			} else if(!strncmp(buffer, "Sbjct", 5)) {
				if(state == between) {
					state = subjectLine;
					alline(buffer, subjectl);
					if(strlen(subjectl) != strlen(queryl))
						throw file_parse_exception(lineNumber);
					get_match(match_mask, queryl, subjectl, record.hit, record.len, record.rid, current_len, record.ungapped_len, q, s, record.stop);
				} else
					throw file_parse_exception(lineNumber);
			} else if(state == queryLine && buffer[0] == ' ') {
				state = between;
			} else if((state == subjectLine || state == separator) && buffer[0] == 10) {
				//printf("%i\n",(int)buffer[0]);
				if(state == subjectLine)
					state = separator;
				else if(state == separator)
					state = end;
				else
					throw file_parse_exception(lineNumber);
			} else {
				if(state != begin && state != queryStart && state != haveid && state != subjectStart)
					throw file_parse_exception(lineNumber);
			}

		}

		if(feof(file))
			return false;
		if(state != end)
			throw file_parse_exception(lineNumber);
		++currentQueryCount;
		++matchCount;
		return true;

	}

	bool get(blast_match &match, const blast_tab_format&)
	{
		char query[nameBufferSize], s2[nameBufferSize];
		char query2[64];
		double f1;
		unsigned long long i1,i2,i3,i4,i5,i6,i7;

		if(readLine(buffer)) {
			while(buffer[0] == '#')
				readLine(buffer);
			if(sscanf(buffer, "%s%s%lf%llu%llu%llu%llu%llu%llu%llu%lf%lf", query, s2, &f1, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &match.expect, &match.bitscore) == 12
				|| sscanf(buffer, "%s%s%s%lf%llu%llu%llu%llu%llu%llu%llu%lf%lf", query, query2, s2, &f1, &i1, &i2, &i3, &i4, &i5, &i6, &i7, &match.expect, &match.bitscore) == 13) {
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
				printf("%s\n", buffer);
				throw file_parse_exception(lineNumber);
			}
		} else {
			if(!text_file::at_end())
				throw file_parse_exception(lineNumber);
			match.set_empty();
			return false;
		}
	}

	template<typename _format>
	void get_read(mcont &v, const _format& format)
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
				return;
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
		return;
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
	char buffer[readBufferSize], currentQuery[nameBufferSize], currentSubject[nameBufferSize];
	blast_match save;

};

#endif /* MATCH_FILE_H_ */
