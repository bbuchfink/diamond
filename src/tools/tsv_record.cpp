#include <stdexcept>
#include "tsv_record.h"

using namespace std;

TextInputFile& operator>>(TextInputFile &file, TSVRecord &record) {
	constexpr size_t BUF_SIZE = 512;
	char query_buf[BUF_SIZE], subject_buf[BUF_SIZE];
	record.qseqid.clear();
	file.getline();
	if (!file.eof()) {
		if (sscanf(file.line.c_str(), "%s%s%lf%zu%zu%zu%zu%zu%zu%zu%lf%lf", query_buf, subject_buf, &record.pident, &record.length, &record.mismatch,
			&record.gapopen, &record.qstart, &record.qend, &record.sstart, &record.send, &record.evalue, &record.bitscore) != 12)
			throw runtime_error("Blast TSV parse error.");
		record.qseqid = query_buf;
		record.sseqid = subject_buf;
	}
	return file;
}

std::ostream& operator<<(std::ostream &str, const TSVRecord &record) {
	str << record.qseqid << '\t'
		<< record.sseqid << '\t'
		<< record.pident << '\t'
		<< record.length << '\t'
		<< record.mismatch << '\t'
		<< record.gapopen << '\t'
		<< record.qstart << '\t'
		<< record.qend << '\t'
		<< record.sstart << '\t'
		<< record.send << '\t'
		<< record.evalue << '\t'
		<< record.bitscore;
	return str;
}