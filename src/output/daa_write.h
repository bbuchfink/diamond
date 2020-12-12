/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef DAA_WRITE_H_
#define DAA_WRITE_H_

#include <limits>
#include <stdint.h>
#include "output.h"
#include "daa_file.h"
#include "../stats/score_matrix.h"
#include "../data/reference.h"
#include "../align/align.h"
#include "../basic/packed_sequence.h"
#include "../data/ref_dictionary.h"
#include "../util/io/consumer.h"

inline void init_daa(OutputFile &f)
{
	DAA_header1 h1;
	f.write(&h1, 1);	
	DAA_header2 h2_;
	f.write(&h2_, 1);
}

inline size_t write_daa_query_record(TextBuffer &buf, const char *query_name, const sequence &query)
{
	size_t seek_pos = buf.size();
	buf.write((uint32_t)0);
	uint32_t l = (uint32_t)query.length();
	buf.write(l);
	buf.write_c_str(query_name, find_first_of(query_name, Const::id_delimiters));
	Packed_sequence s(query, align_mode.input_sequence_type);
	uint8_t flags = s.has_n() ? 1 : 0;
	buf.write(flags);
	buf << s.data();
	return seek_pos;
}

inline void finish_daa_query_record(TextBuffer &buf, size_t seek_pos)
{
	*(uint32_t*)(&buf[seek_pos]) = (uint32_t)(buf.size() - seek_pos - sizeof(uint32_t));
}

inline void write_daa_record(TextBuffer &buf, const IntermediateRecord &r)
{
	buf.write(r.subject_dict_id).write(r.flag);
	buf.write_packed(r.score);
	buf.write_packed(r.query_begin);
	buf.write_packed(r.subject_begin);
	buf << r.transcript.data();
}

inline void write_daa_record(TextBuffer &buf, const Hsp &match, size_t subject_id)
{
	buf.write(config.command == Config::view ? (uint32_t)subject_id : ReferenceDictionary::get().get(current_ref_block, subject_id));
	buf.write(get_segment_flag(match));
	buf.write_packed(match.score);
	buf.write_packed(match.oriented_range().begin_);
	buf.write_packed(match.subject_range.begin_);
	buf << match.transcript.data();
}

inline void finish_daa(OutputFile &f, const DatabaseFile &db)
{
	DAA_header2 h2_(db.ref_header.sequences,
		config.db_size,
		score_matrix.gap_open(),
		score_matrix.gap_extend(),
		config.reward,
		config.penalty,
		score_matrix.k(),
		score_matrix.lambda(),
		config.max_evalue,
		to_lower_case(config.matrix),
		align_mode.mode);

	h2_.block_type[0] = DAA_header2::alignments;
	h2_.block_type[1] = DAA_header2::ref_names;
	h2_.block_type[2] = DAA_header2::ref_lengths;

	const ReferenceDictionary &dict = ReferenceDictionary::get();

	uint32_t size = 0;
	f.write(&size, 1);
	h2_.block_size[0] = f.tell() - sizeof(DAA_header1) - sizeof(DAA_header2);
	h2_.db_seqs_used = dict.seqs();
	h2_.query_records = statistics.get(Statistics::ALIGNED);

	size_t s = 0;
	for (PtrVector<string>::const_iterator i = dict.name_.begin(); i != dict.name_.end(); ++i) {
		//std::cout << **i << endl;
		f << **i;
		s += (*i)->length() + 1;
	}
	h2_.block_size[1] = s;

	f.write(dict.len_.data(), dict.len_.size());
	h2_.block_size[2] = dict.len_.size() * sizeof(uint32_t);

	f.seek(sizeof(DAA_header1));
	f.write(&h2_, 1);
}

inline void finish_daa(OutputFile &f, DAA_file &daa_in) {
	DAA_header2 h2_(daa_in.db_seqs(),
		daa_in.db_letters(),
		daa_in.gap_open_penalty(),
		daa_in.gap_extension_penalty(),
		daa_in.match_reward(),
		daa_in.mismatch_penalty(),
		daa_in.kappa(),
		daa_in.lambda(),
		daa_in.evalue(),
		daa_in.score_matrix(),
		daa_in.mode());

	h2_.block_type[0] = DAA_header2::alignments;
	h2_.block_type[1] = DAA_header2::ref_names;
	h2_.block_type[2] = DAA_header2::ref_lengths;

	uint32_t size = 0;
	f.write(&size, 1);
	h2_.block_size[0] = f.tell() - sizeof(DAA_header1) - sizeof(DAA_header2);
	h2_.db_seqs_used = daa_in.db_seqs_used();
	h2_.query_records = daa_in.query_records();

	for (size_t i = 0; i < daa_in.db_seqs_used(); ++i)
		f << daa_in.ref_name(i);
	h2_.block_size[1] = daa_in.block_size(1);

	f.write(daa_in.ref_len().data(), daa_in.ref_len().size());
	h2_.block_size[2] = daa_in.block_size(2);

	f.seek(sizeof(DAA_header1));
	f.write(&h2_, 1);
}

#endif /* DAA_WRITE_H_ */
