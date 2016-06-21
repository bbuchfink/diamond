#ifndef OUTPUT_H_
#define OUTPUT_H_

#include "../util/binary_file.h"
#include "../basic/packed_transcript.h"

struct Intermediate_record
{
	void read(Buffered_file &f)
	{
		f.read(query_id);
		f.read(subject_id);
		f.read(flag);
		f.read_packed(flag & 3, score);
		f.read_packed((flag >> 2) & 3, query_begin);
		f.read_packed((flag >> 4) & 3, subject_begin);
		transcript.read(f);
	}
	uint32_t query_id, subject_id, score, query_begin, subject_begin;
	uint8_t flag;
	Packed_transcript transcript;
};

#endif