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

#include <fstream>
#include "fasta_file.h"
#include "../../util/system/system.h"
#include "../../util/seq_file_format.h"
#include "../../basic/config.h"
#include "../../util/sequence/sequence.h"
#include "../../util/string/tokenizer.h"

using std::vector;
using std::string;
using std::pair;
using std::runtime_error;
using std::tie;
using namespace Util::Tsv;

FastaFile::FastaFile(const std::vector<std::string>& file_name, Metadata metadata, Flags flags, const ValueTraits& value_traits):
	SequenceFile(SequenceFile::Type::FASTA, Alphabet::STD, flags, FormatFlags::DICT_LENGTHS | FormatFlags::DICT_SEQIDS, value_traits),
	oid_(0)
{
	if (file_name.size() > 2)
		throw OperationNotSupported();
	for (const auto& f : file_name)
		file_.emplace_back(f);
	file_ptr_ = file_.begin();
	format_ = guess_format(file_.front());
	if (!flag_any(flags, Flags::NEED_LETTER_COUNT))
		return;
#if EXTRA
	/*if (exists(file_name.front() + ".fai")) {
		task_timer timer("Reading faidx file", 3);
		int64_t seqs = 0, letters = 0;
		string acc;
		for (const auto& f : file_name) {
			pair<int64_t, int64_t> r = read_fai_file(f + ".fai", seqs, letters);
			seqs += r.first;
			letters += r.second;
		}
		seqs_ = seqs;
		letters_ = letters;
	}
	else {*/
#endif
		task_timer timer("Reading fasta file", 3);
		tie(seqs_, letters_) = init_read();
		set_seqinfo_ptr(0);
#ifdef EXTRA
	//}
#endif
}

FastaFile::FastaFile(const string& file_name, bool overwrite, const WriteAccess&):
	SequenceFile(SequenceFile::Type::FASTA, Alphabet::STD, Flags::NONE, FormatFlags::DICT_LENGTHS | FormatFlags::DICT_SEQIDS, amino_acid_traits),
	out_file_(file_name.empty() ? new TempFile : new OutputFile(file_name, Compressor::NONE, overwrite ? "w+b" : "r+b")),
	format_(new FASTA_format()),
	oid_(0),
	seqs_(0),
	letters_(0)
{
	file_.emplace_back(*out_file_);
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
		if (++file_ptr_ != file_.end() && format_->get_seq(id, seq, *file_ptr_, value_traits_, nullptr))
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
			if (f.temp_file)
				f.close_and_delete();
			else
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
	const bool r = format_->get_seq(id, seq, *(file_ptr_++), this->value_traits_, quals);
	if (file_ptr_ == file_.end())
		file_ptr_ = file_.begin();
	return r;
}

void FastaFile::create_partition_balanced(size_t max_letters) {
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

size_t FastaFile::letters() const {
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

void FastaFile::close_weakly()
{
}

void FastaFile::reopen()
{
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
	return file_.front().file_name;
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
	throw OperationNotSupported();
}

void FastaFile::init_random_access(const size_t query_block, const size_t ref_blocks, bool dictionary)
{
	if (dictionary)
		load_dictionary(query_block, ref_blocks);
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
	Util::Seq::format(seq, id.c_str(), nullptr,  *out_file_, "fasta", value_traits_);
	++seqs_;
	letters_ += seq.length();
}

void FastaFile::prep_db(const string& path) {
	task_timer timer("Indexing FASTA file");
	TextInputFile f(path);
	FASTA_format fmt;
	vector<Letter> seq;
	string id;
	std::ofstream out(path + ".fai");
	int64_t seqs = 0, letters = 0;
	while (fmt.get_seq(id, seq, f, amino_acid_traits)) {
		out << Util::Seq::seqid(id.c_str(), false) << '\t' << seq.size() << std::endl;
		++seqs;
		letters += seq.size();
	}
	f.close();
	out.close();
	timer.finish();
	message_stream << "Processed " << seqs << " sequences, " << letters << " total letters." << std::endl;
}

int64_t FastaFile::line_count() const {
	return file_.front().line_count;
}

std::pair<int64_t, int64_t> FastaFile::init_read() {
	vector<Letter> seq;
	string id;
	int64_t seqs = 0, letters = 0;
	while (read_seq(seq, id)) {
		if (flag_any(flags_, Flags::ACC_TO_OID_MAPPING | Flags::OID_TO_ACC_MAPPING))
			add_seqid_mapping(id, seqs);
		++seqs;
		letters += seq.size();
	}
	return { seqs, letters };
}