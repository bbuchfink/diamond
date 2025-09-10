/****
DIAMOND protein aligner
Copyright (C) 2021 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#include "fasta_file.h"
#include "util/sequence/sequence.h"

using std::vector;
using std::string;
using std::pair;
using std::runtime_error;
using std::unique_ptr;
using std::tie;
using std::next;
using namespace Util::Tsv;
using Util::String::TokenizerBase;
using Util::String::FastaTokenizer;
using Util::String::FastqTokenizer;

static const Schema FASTA_SCHEMA { Type::STRING, Type::STRING };
static const Schema FASTQ_SCHEMA{ Type::STRING, Type::STRING, Type::STRING };
//const char* FASTA_SEP = "\n>", * FASTQ_SEP = "\n@";

SeqFileFormat guess_format(TextInputFile& file) {
	string r = file.peek(1);
	if (r.empty())
		throw std::runtime_error("Error detecting input file format. Input file seems to be empty.");
	switch (r.front()) {
	case '>': return SeqFileFormat::FASTA;
	case '@': return SeqFileFormat::FASTQ;
	default: throw std::runtime_error("Error detecting input file format. First line must begin with '>' (FASTA) or '@' (FASTQ).");
	}
}

FastaFile::FastaFile(const std::vector<std::string>& file_name, Metadata metadata, Flags flags, const ValueTraits& value_traits):
	SequenceFile(SequenceFile::Type::FASTA, Alphabet::STD, flags, FormatFlags::DICT_LENGTHS | FormatFlags::DICT_SEQIDS, value_traits),
	oid_(0)
{
	if (file_name.size() > 2)
		throw OperationNotSupported();
	unique_ptr<TextInputFile> input_file(new TextInputFile(file_name.front()));
	format_ = guess_format(*input_file);
	//Util::Tsv::Config config(format_ == SeqFileFormat::FASTA ? FASTA_SEP : FASTQ_SEP, format_ == SeqFileFormat::FASTA ? (TokenizerBase*)(new FastaTokenizer()) : new FastqTokenizer);
	Util::Tsv::Config config(format_ == SeqFileFormat::FASTA ? (TokenizerBase*)(new FastaTokenizer()) : new FastqTokenizer);
	const auto schema = format_ == SeqFileFormat::FASTA ? FASTA_SCHEMA : FASTQ_SCHEMA;
	file_.emplace_back(schema, std::move(input_file), Util::Tsv::Flags(), config);
	if (file_name.size() > 1)
		file_.emplace_back(schema, file_name[1], Util::Tsv::Flags(), config);
	file_ptr_ = file_.begin();
	
	if (!flag_any(flags, Flags::NEED_LETTER_COUNT))
		return;
	tie(seqs_, letters_) = init_read();
	set_seqinfo_ptr(0);
}

FastaFile::FastaFile(const string& file_name, bool overwrite, const WriteAccess&, Flags flags, const ValueTraits& value_traits):
	SequenceFile(SequenceFile::Type::FASTA, Alphabet::STD, flags, FormatFlags::DICT_LENGTHS | FormatFlags::DICT_SEQIDS, value_traits),
	out_file_(file_name.empty() ? new TempFile : new OutputFile(file_name, Compressor::NONE, overwrite ? "w+b" : "r+b")),
	format_(SeqFileFormat::FASTA),
	oid_(0),
	seqs_(0),
	letters_(0)
{
	unique_ptr<TextInputFile> in(new TextInputFile(*out_file_));
	//Util::Tsv::Config config(FASTA_SEP, new FastaTokenizer());
	Util::Tsv::Config config(new FastaTokenizer());
	file_.emplace_back(FASTA_SCHEMA, std::move(in), Util::Tsv::Flags(), config);
	//file_.emplace_back(*out_file_);
	file_ptr_ = file_.begin();
	if (!overwrite) {
		vector<Letter> seq;
		string id;
		while (read_seq(seq, id, nullptr)) {
			++seqs_;
			letters_ += seq.size();
		}
	}
}

int64_t FastaFile::file_count() const {
	return file_.size();
}

bool FastaFile::files_synced() {
	if (file_ptr_ != file_.begin())
		return false;
	if (file_ptr_->eof()) {
		string id;
		vector<Letter> seq;
		if (next(file_ptr_) != file_.end() && !next(file_ptr_)->read_record().empty())
			return false;
	}
	return true;
}

SequenceFile::SeqInfo FastaFile::read_seqinfo() {
	throw OperationNotSupported();
}

void FastaFile::putback_seqinfo() {
	throw OperationNotSupported();
}

void FastaFile::close() {
	if (out_file_)
		file_.front().close();
	else
		for (auto& f : file_)
			f.close();
	
}

void FastaFile::set_seqinfo_ptr(OId i) {
	if (out_file_)
		out_file_->rewind();
	for (auto& f : file_)
		f.rewind();
	oid_ = 0;
	vector<Letter> seq;
	string id;
	while (oid_ != i) {
		read_seq(seq, id, nullptr);
	}
}

OId FastaFile::tell_seq() const {
	return oid_;
}

bool FastaFile::eof() const {
	return file_.back().eof();
}

void FastaFile::init_seq_access() {
	set_seqinfo_ptr(0);
}

bool FastaFile::read_seq(vector<Letter>& seq, string &id, std::vector<char>* quals)
{
	oid_++;
	Table t = file_ptr_->read_record();
	if (!t.empty()) {
		id = t.front().get<string>(0);
		Util::Seq::from_string(t.front().get<string>(1), seq, value_traits_, 0);
		if (format_ == SeqFileFormat::FASTQ && quals) {
			const string q = t.front().get<string>(2);
			quals->assign(q.begin(), q.end());
		}
		if (++file_ptr_ == file_.end())
			file_ptr_ = file_.begin();
	}
	return !t.empty();
}

void FastaFile::create_partition_balanced(int64_t max_letters) {
	throw OperationNotSupported();
}

void FastaFile::save_partition(const string & partition_file_name, const string & annotation) {
	throw OperationNotSupported();
}

int FastaFile::get_n_partition_chunks() {
	throw OperationNotSupported();
}

void FastaFile::init_seqinfo_access() {
	throw OperationNotSupported();
}

void FastaFile::seek_chunk(const Chunk& chunk) {
	throw OperationNotSupported();
}

size_t FastaFile::id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) {
	throw OperationNotSupported();
}

void FastaFile::seek_offset(size_t p) {
}

void FastaFile::read_seq_data(Letter* dst, size_t len, size_t& pos, bool seek) {
	throw OperationNotSupported();
}

void FastaFile::read_id_data(const int64_t oid, char* dst, size_t len) {
	throw OperationNotSupported();
}

void FastaFile::skip_id_data() {
	throw OperationNotSupported();
}

int64_t FastaFile::sequence_count() const {
	return seqs_;
}

int64_t FastaFile::letters() const {
	return letters_;
}

int FastaFile::db_version() const {
	throw OperationNotSupported();
}

int FastaFile::program_build_version() const {
	throw OperationNotSupported();
}

SequenceFile::Metadata FastaFile::metadata() const {
	return Metadata();
}

int FastaFile::build_version() {
	throw OperationNotSupported();
}

FastaFile::~FastaFile()
{
	close();
}

BitVector* FastaFile::filter_by_accession(const std::string& file_name)
{
	throw std::runtime_error("The FASTA database format does not support filtering by accession.");
}

const BitVector* FastaFile::builtin_filter()
{
	return nullptr;
}

std::string FastaFile::file_name()
{
	return file_.front().file_name();
}

int64_t FastaFile::sparse_sequence_count() const
{
	throw OperationNotSupported();
}

std::vector<TaxId> FastaFile::taxids(size_t oid) const
{
	throw OperationNotSupported();
}

void FastaFile::seq_data(size_t oid, std::vector<Letter>& dst) const
{
	throw OperationNotSupported();
}

size_t FastaFile::seq_length(size_t oid) const
{
	if (oid < seq_length_.size())
		return seq_length_[oid];
	throw OperationNotSupported();
}

void FastaFile::end_random_access(bool dictionary)
{
	if (!dictionary)
		return;
	free_dictionary();
}

void FastaFile::init_write() {
	out_file_->seek(0, SEEK_END);
}

void FastaFile::write_seq(const Sequence& seq, const std::string& id) {
	static TextBuffer buf;
	Util::Seq::format(seq, id.c_str(), nullptr,  buf, "fasta", value_traits_);
	out_file_->write(buf.data(), buf.size());
	buf.clear();
	++seqs_;
	letters_ += seq.length();
	if (flag_any(flags_, Flags::NEED_LENGTH_LOOKUP))
		seq_length_.push_back(seq.length());
}

std::pair<int64_t, int64_t> FastaFile::init_read() {
	vector<Letter> seq;
	string id;
	int64_t seqs = 0, letters = 0;
	while (read_seq(seq, id)) {
		if (flag_any(flags_, Flags::ACC_TO_OID_MAPPING | Flags::OID_TO_ACC_MAPPING))
			add_seqid_mapping(id, seqs);
		if (flag_any(flags_, Flags::NEED_LENGTH_LOOKUP))
			seq_length_.push_back((Loc)seq.size());
		++seqs;
		letters += seq.size();
	}
	return { seqs, letters };
}