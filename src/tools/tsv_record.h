#ifndef TSV_RECORD_H_
#define TSV_RECORD_H_

#include <string>
#include <ostream>
#include "../util/io/text_input_file.h"

struct TSVRecord  {
	std::string qseqid, sseqid;
	double pident, evalue, bitscore;
	size_t length, mismatch, gapopen, qstart, qend, sstart, send;
	
	bool blank() const {
		return qseqid.empty();
	}

	friend TextInputFile& operator>>(TextInputFile &file, TSVRecord &record);
	friend std::ostream& operator<<(std::ostream &str, const TSVRecord &record);
};

#endif