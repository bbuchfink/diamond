/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#include "fasta_file.h"
#include "util/sequence/sequence.h"
#include "util/system/system.h"

using std::vector;
using std::string;
using std::pair;
using std::runtime_error;
using std::unique_ptr;
using std::tie;
using std::ifstream;
using std::ofstream;
using std::next;
using std::unordered_map;
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

FastaFile::FastaFile(const vector<string>& file_name, Flags flags, const ValueTraits& value_traits, const string& index_file):
	SequenceFile(SequenceFile::Type::FASTA, flags, FormatFlags::DICT_LENGTHS | FormatFlags::DICT_SEQIDS, value_traits),
	index_file_(index_file),
	oid_(0),
	raw_chunk_no_(0),
	seqs_(0),
	letters_(0)
{
	if (file_name.size() > 2)
		throw OperationNotSupported();
	if (bool(flags & (Flags::TAXON_MAPPING | Flags::TAXON_NODES | Flags::TAXON_RANKS | Flags::TAXON_SCIENTIFIC_NAMES)))
		throw runtime_error("Fasta database format does not support taxonomic features.");
	unique_ptr<TextInputFile> input_file(new TextInputFile(file_name.front()));
	format_ = guess_format(*input_file);
	//Util::Tsv::Config config(format_ == SeqFileFormat::FASTA ? FASTA_SEP : FASTQ_SEP, format_ == SeqFileFormat::FASTA ? (TokenizerBase*)(new FastaTokenizer()) : new FastqTokenizer);
	Util::Tsv::Config cfg(format_ == SeqFileFormat::FASTA ? (TokenizerBase*)(new FastaTokenizer()) : new FastqTokenizer);
	const auto schema = format_ == SeqFileFormat::FASTA ? FASTA_SCHEMA : FASTQ_SCHEMA;
	file_.emplace_back(schema, std::move(input_file), Util::Tsv::Flags(), cfg);
	if (file_name.size() > 1)
		file_.emplace_back(schema, file_name[1], Util::Tsv::Flags(), cfg);
	file_ptr_ = file_.begin();
	if (!config.fasta_index_file.empty()) {
		ifstream in(config.fasta_index_file, std::ios::binary);
		uint64_t pos;
		while (in >> pos)
			index_.push_back(pos);
	}
	if (!flag_any(flags, Flags::NEED_LETTER_COUNT))
		return;
	tie(seqs_, letters_) = init_read();
	set_seqinfo_ptr(0);
}

FastaFile::FastaFile(const string& file_name, bool overwrite, const WriteAccess&, Flags flags, const ValueTraits& value_traits):
	SequenceFile(SequenceFile::Type::FASTA, flags, FormatFlags::DICT_LENGTHS | FormatFlags::DICT_SEQIDS, value_traits),
	out_file_(file_name.empty() ? new TempFile : new OutputFile(file_name, Compressor::NONE, overwrite ? "w+b" : "r+b")),
	format_(SeqFileFormat::FASTA),
	oid_(0),
	raw_chunk_no_(0),
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
	if (i == oid_)
		return;
	if (i == 0) {
		raw_chunk_no_ = 0;
		oid_ = 0;
		if (out_file_)
			out_file_->rewind();
		for (auto& f : file_)
			f.rewind();
		return;
	}
	else if (i == seqs_) {
		oid_ = seqs_;
		for (auto& f : file_)
			f.seek(f.size());
		return;
	}
	if (index_.empty()) {
		throw runtime_error("FastaFile::set_seqinfo_ptr: no index available");
	}
	else {
		if (file_.size() != 1 || i >= index_.size())
			throw runtime_error("FastaFile::set_seqinfo_ptr");
		file_.front().seek(index_[i]);
		oid_ = i;
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
			const string q = Util::Seq::remove_newlines(t.front().get<string>(2));
			quals->assign(q.begin(), q.end());
		}
		if (++file_ptr_ == file_.end())
			file_ptr_ = file_.begin();
	}
	return !t.empty();
}

namespace {

struct FastaRawChunk : public RawChunk {
	FastaRawChunk(const ValueTraits& value_traits) :
		value_traits_(value_traits),
		bytes_(0)
	{}

	virtual bool empty() const override {
		return data_.empty();
	}

	virtual OId begin() const noexcept override {
		return 0;
	}

	virtual OId end() const noexcept override {
		return 0;
	}

	virtual DecodedPackage* decode(SequenceFile::Flags flags, const BitVector* filter, unordered_map<string, bool>* accs) const override {
		assert(filter == nullptr && accs == nullptr);
		DecodedPackage* pkg = new DecodedPackage();
		pkg->no = no;
		const bool titles = flag_any(flags, SequenceFile::Flags::TITLES),
			seqs = flag_any(flags, SequenceFile::Flags::SEQS);
		FastaTokenizer tokenizer;
		vector<Letter> seq;
		const char* ptr = data_.data(), *end = ptr + data_.size();
		while (ptr < end) {
			const char* record_end = find_record_end(ptr, end);
			Table table(FASTA_SCHEMA);
			table.push_back(ptr, record_end, &tokenizer);
			if (seqs) {
				Util::Seq::from_string(table.front().get<string>(1), seq, value_traits_, 0);
				pkg->seqs.push_back(seq.begin(), seq.end());
			}
			if (titles) {
				const string id = table.front().get<string>(0);
				pkg->ids.push_back(id.begin(), id.end());				
			}
			++pkg->seq_count;
			ptr = record_end == end ? end : record_end + 1;
		}
		return pkg;
	}

	virtual size_t letters() const noexcept override {
		return 0;
	}

	virtual size_t bytes() const noexcept override {
		return bytes_;
	}

	static const char* find_record_end(const char* begin, const char* end) {
		static const char* delimiter = "\n>";
		return std::search(begin + 1, end, delimiter, delimiter + 2);
	}

	const ValueTraits& value_traits_;
	string data_;
	size_t bytes_;
};

}

RawChunk* FastaFile::raw_chunk(size_t bytes, Flags flags) {
	if (format_ != SeqFileFormat::FASTA)
		throw OperationNotSupported();
	FastaRawChunk* chunk = new FastaRawChunk(value_traits_);
	chunk->no = raw_chunk_no_++;
	if (bytes == 0)
		return chunk;
	chunk->bytes_ = file_ptr_->read_raw(chunk->data_, bytes);
	if (chunk->data_.empty())
		return chunk;
	bool complete = chunk->bytes_ < bytes;
	if (!complete && chunk->data_.back() == '\n') {
		const string next = file_ptr_->peek(1);
		complete = next.empty() || next.front() == '>';
	}
	if (!complete) {
		const pair<bool, size_t> r = file_ptr_->read_to_fasta_record_end(chunk->data_);
		chunk->bytes_ += r.second;
	}
	(void)flags;
	return chunk;
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

void FastaFile::read_id_data(const int64_t oid, char* dst, size_t len, bool all, bool full_titles) {
	throw OperationNotSupported();
}

void FastaFile::skip_id_data() {
	throw OperationNotSupported();
}

uint64_t FastaFile::sequence_count() const {
	return seqs_;
}

uint64_t FastaFile::letters() const {
	return letters_;
}

int FastaFile::db_version() const {
	throw OperationNotSupported();
}

int FastaFile::program_build_version() const {
	throw OperationNotSupported();
}

int FastaFile::build_version() {
	throw OperationNotSupported();
}

FastaFile::~FastaFile()
{
	close();
}

DbFilter* FastaFile::filter_by_accession(const std::string& file_name)
{
	throw std::runtime_error("The FASTA database format does not support filtering by accession.");
}

std::string FastaFile::file_name()
{
	return file_.front().file_name();
}

std::vector<TaxId> FastaFile::taxids(size_t oid) const
{
	throw OperationNotSupported();
}

void FastaFile::seq_data(size_t oid, std::vector<Letter>& dst)
{
	throw OperationNotSupported();
}

Loc FastaFile::seq_length(size_t oid)
{
	if (oid < seq_length_.size())
		return seq_length_[oid];
	throw runtime_error("Missing FASTA index");
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

pair<int64_t, int64_t> FastaFile::init_read() {
	int64_t seqs = 0, letters = 0;
	const bool count_only = !flag_any(flags_, Flags::ACC_TO_OID_MAPPING | Flags::OID_TO_ACC_MAPPING | Flags::NEED_LENGTH_LOOKUP);
	if (index_file_.empty() || !exists(index_file_)) {
		const uint64_t limit = std::min<uint64_t>(Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT)), SequenceFile::DEFAULT_LOAD_SIZE);
		uint64_t bytes = 0, ms = 0;
		OId oid = 0;
		const SequenceFile::Flags flags = flags_;
		flags_ |= Flags::SEQS;
		if(flag_any(flags_, Flags::ACC_TO_OID_MAPPING | Flags::OID_TO_ACC_MAPPING))
			flags_ |= Flags::TITLES;
		for (;;) {
			TaskTimer timer;
			Block* b = load_seqs(limit);
			ms += timer.microseconds();
			bytes += b->raw_bytes();
			if (b->empty()) {
				delete b;
				break;
			}
			const BlockId n = b->seqs().size();
			letters += b->seqs().letters();
			seqs += n;
			if (count_only) {
				delete b;
				continue;
			}
			seq_length_.reserve(seqs);
			for (BlockId i = 0; i < n; ++i) {
				if (flag_any(flags_, Flags::ACC_TO_OID_MAPPING | Flags::OID_TO_ACC_MAPPING))
					add_seqid_mapping(b->ids()[i], oid);
				if (flag_any(flags_, Flags::NEED_LENGTH_LOOKUP))
					seq_length_.push_back(b->seqs().length(i));
				++oid;
			}
			delete b;
		}
		std::ostringstream ss;
		ss << "Loaded " << bytes << " bytes from disk at " << ((double)bytes / MEGABYTES / ms * 1e6) << " MB/s" << std::endl;
		open_stats_ = ss.str();
		if (flag_any(flags_, Flags::NEED_LENGTH_LOOKUP) && !index_file_.empty()) {
			ofstream out(index_file_, std::ios::binary);
			out.write((const char*)seq_length_.data(), seq_length_.size() * sizeof(Loc));
			if (!out)
				throw runtime_error("Error writing file: " + index_file_);
		}
		flags_ = flags;
	}
	else {
		ifstream in(index_file_, std::ios::binary);
		if (!in)
			throw runtime_error("Error opening file: " + index_file_);
		in.seekg(0, std::ios::end);
		seqs = in.tellg() / sizeof(Loc);
		seq_length_.resize(seqs);
		in.seekg(0, std::ios::beg);
		in.read((char*)seq_length_.data(), seqs * sizeof(Loc));
		for (size_t i = 0; i < seqs; ++i)
			letters += seq_length_[i];
	}
	return { seqs, letters };
}