#include "daa_write.h"
#include "../util/util.h"
#include "../util/sequence/sequence.h"
#include "../../basic/statistics.h"

using std::string;
using std::vector;

void init_daa(OutputFile& f)
{
	DAA_header1 h1;
	f.write(&h1, 1);
	DAA_header2 h2_;
	f.write(&h2_, 1);
}

size_t write_daa_query_record(TextBuffer& buf, const char* query_name, const Sequence& query)
{
	size_t seek_pos = buf.size();
	buf.write((uint32_t)0);
	uint32_t l = (uint32_t)query.length();
	buf.write(l);
	buf.write_c_str(query_name, find_first_of(query_name, Util::Seq::id_delimiters));
	Packed_sequence s(query, align_mode.input_sequence_type);
	uint8_t flags = s.has_n() ? 1 : 0;
	buf.write(flags);
	buf << s.data();
	return seek_pos;
}

void finish_daa_query_record(TextBuffer& buf, size_t seek_pos)
{
	*(uint32_t*)(&buf[seek_pos]) = (uint32_t)(buf.size() - seek_pos - sizeof(uint32_t));
}

void write_daa_record(TextBuffer& buf, const IntermediateRecord& r)
{
	buf.write((uint32_t)r.target_dict_id).write(r.flag);
	buf.write_packed(r.score);
	buf.write_packed(r.query_begin);
	buf.write_packed(r.subject_begin);
	buf << r.transcript.data();
}

void write_daa_record(TextBuffer& buf, const Hsp& match, uint32_t subject_id)
{
	buf.write(subject_id);
	buf.write(get_segment_flag(match));
	buf.write_packed(match.score);
	buf.write_packed(match.oriented_range().begin_);
	buf.write_packed(match.subject_range.begin_);
	buf << match.transcript.data();
}

void finish_daa(OutputFile& f, const SequenceFile& db)
{
	DAA_header2 h2_(db.sequence_count(),
		score_matrix.db_letters(),
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

	uint32_t size = 0;
	f.write(&size, 1);
	h2_.block_size[0] = f.tell() - sizeof(DAA_header1) - sizeof(DAA_header2);
	const size_t n = db.dict_size();
	h2_.db_seqs_used = n;
	h2_.query_records = statistics.get(Statistics::ALIGNED);

	size_t s = 0;
	for (size_t i = 0; i < n; ++i) {
		const string title = db.dict_title(i, 0);
		f << title;
		s += title.length() + 1;
	}
	h2_.block_size[1] = s;

	for (size_t i = 0; i < n; ++i)
		f << (uint32_t)db.dict_len(i, 0);
	h2_.block_size[2] = n * sizeof(uint32_t);

	f.seek(sizeof(DAA_header1));
	f.write(&h2_, 1);
}

DAA_header2::DAA_header2(const DAA_file& daa_in) :
	DAA_header2(daa_in.db_seqs(),
		daa_in.db_letters(),
		daa_in.gap_open_penalty(),
		daa_in.gap_extension_penalty(),
		daa_in.match_reward(),
		daa_in.mismatch_penalty(),
		daa_in.kappa(),
		daa_in.lambda(),
		daa_in.evalue(),
		daa_in.score_matrix(),
		daa_in.mode())
{
	block_type[0] = DAA_header2::alignments;
	block_type[1] = DAA_header2::ref_names;
	block_type[2] = DAA_header2::ref_lengths;
}

static void terminate_aln_block(OutputFile& f, DAA_header2& h2) {
	uint32_t size = 0;
	f.write(&size, 1);
	h2.block_size[0] = f.tell() - sizeof(DAA_header1) - sizeof(DAA_header2);
}

static void write_header2(OutputFile& f, const DAA_header2& h2) {
	f.seek(sizeof(DAA_header1));
	f.write(&h2, 1);
}

void finish_daa(OutputFile& f, DAA_file& daa_in) {
	DAA_header2 h2_(daa_in);

	terminate_aln_block(f, h2_);
	
	h2_.db_seqs_used = daa_in.db_seqs_used();
	h2_.query_records = daa_in.query_records();

	for (size_t i = 0; i < daa_in.db_seqs_used(); ++i)
		f << daa_in.ref_name(i);
	h2_.block_size[1] = daa_in.block_size(1);

	f.write(daa_in.ref_len().data(), daa_in.ref_len().size());
	h2_.block_size[2] = daa_in.block_size(2);

	write_header2(f, h2_);
}

void finish_daa(OutputFile& f, DAA_file& daa_in, const StringSet& seq_ids, const vector<uint32_t>& seq_lens, int64_t query_count) {
	DAA_header2 h2_(daa_in);
	terminate_aln_block(f, h2_);
	h2_.db_seqs_used = seq_ids.size();
	h2_.query_records = (uint64_t)query_count;
	int64_t s = 0;
	for (int64_t i = 0; i < seq_ids.size(); ++i) {
		f << seq_ids[i];
		s += seq_ids.length(i) + 1;
	}
	h2_.block_size[1] = (uint64_t)s;
	f.write(seq_lens.data(), seq_lens.size());
	h2_.block_size[2] = seq_lens.size() * sizeof(uint32_t);
	write_header2(f, h2_);
}
