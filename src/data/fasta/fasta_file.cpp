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

#include <fstream>
#include <thread>
#include "fasta_file.h"
#include "util/sequence/sequence.h"
#include "util/system/system.h"
#include "util/log_stream.h"
#include "parser.h"

using std::vector;
using std::string;
using std::pair;
using std::runtime_error;
using std::unique_ptr;
using std::tie;
using std::ifstream;
using std::ofstream;
using std::next;
using std::endl;
using std::unordered_map;
using Util::String::FastqTokenizer;

static constexpr int64_t CHECK_FOR_DNA_COUNT = 10;
//const char* FASTA_SEP = "\n>", * FASTQ_SEP = "\n@";

static SeqFileFormat guess_format(File& file) {
	string r = file.peek(1);
	if (r.empty())
		throw EmptyFileError("Error detecting input file format. Input file seems to be empty.");
	switch (r.front()) {
	case '>': return SeqFileFormat::FASTA;
	case '@': return SeqFileFormat::FASTQ;
	default: throw FormatDetectionError("Error detecting input file format (when treating as text file). First line must begin with '>' (FASTA) or '@' (FASTQ).");
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
	
	file_.emplace_back(file_name.front(), "rb", File::Flags::TREAT_BLANK_AS_STDIN | File::Flags::DETECT_COMPRESSION);
	format_ = guess_format(file_.back());
	if (file_name.size() > 1)
		file_.emplace_back(file_name[1], "rb", File::Flags::TREAT_BLANK_AS_STDIN | File::Flags::DETECT_COMPRESSION);
	file_ptr_ = file_.begin();
	if (!config.fasta_index_file.empty()) {
		ifstream in(config.fasta_index_file, std::ios::binary);
		uint64_t pos;
		while (in >> pos)
			index_.push_back(pos);
	}
	size_t size = 0;
	for (const string& name : file_name) {
		if (name.empty() || name == "-")
			continue;
		size += ::file_size(name.c_str());
	}
	disk_size_ = size;
	std::ostringstream ss;
	ss << "FASTA file size: " << disk_size_ << " bytes" << endl;
	open_stats_ += ss.str();
	if (!flag_any(flags, Flags::NEED_LETTER_COUNT))
		return;
	tie(seqs_, letters_) = init_read();
	set_seqinfo_ptr(0);
}

int64_t FastaFile::file_count() const {
	return file_.size();
}

SequenceFile::SeqInfo FastaFile::read_seqinfo() {
	throw OperationNotSupported();
}

void FastaFile::putback_seqinfo() {
	throw OperationNotSupported();
}

void FastaFile::print_info() const {
	*message_stream << "Database: " << config.database << ' ';
	*message_stream << "(type: " << to_string(type());
	if (seqs_ > 0) {
		*message_stream << ", ";
		*message_stream << "sequences: " << sequence_count().value() << ", ";
		*message_stream << "letters: " << letters().value();
	}
	*message_stream << ')' << endl;
}

void FastaFile::close() {
	for (auto& f : file_)
		f.close();	
}

void FastaFile::rewind() {
	raw_chunk_no_ = 0;
	oid_ = 0;
	for (auto& f : file_)
		f.rewind();
	file_ptr_ = file_.begin();
}

void FastaFile::set_seqinfo_ptr(OId i) {
	if (i == oid_)
		return;
	if (i == 0) {
		rewind();
		return;
	}
	else if (i == seqs_) {
		oid_ = seqs_;
		for (auto& f : file_)
			f.seek(0, SEEK_END);
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
	bool ok;
	if (format_ == SeqFileFormat::FASTQ) {		
		FastqTokenizer tok;
		ok = tok.read_record(*file_ptr_, id, seq, quals);
	}
	else {
		ok = Util::String::read_fasta_record(*file_ptr_, id, seq);
	}
	if (!ok)
		return false;
	for (Letter& c : seq)
		c = value_traits_.from_char(c);
	if (++file_ptr_ == file_.end())
		file_ptr_ = file_.begin();
	return true;
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

	virtual DecodedPackage* decode(SequenceFile::Flags flags, const BitVector* filter, unordered_map<string, bool>* accs, SequenceType seq_type) const override {
		static const char* DNA_ERR = "The sequences are expected to be proteins but only contain DNA letters. Use the option --ignore-warnings to proceed.";
		assert(filter == nullptr && accs == nullptr);
		const int frame_mask = ::frame_mask();
		DecodedPackage* pkg = new DecodedPackage();
		pkg->no = no;
		const bool titles = flag_any(flags, SequenceFile::Flags::TITLES),
			seqs = flag_any(flags, SequenceFile::Flags::SEQS);
		vector<Letter> seq;
		string id;
		static constexpr size_t SEQ_LEN_EST = 200;
		const size_t record_est = std::max<size_t>(1, data_.size() / SEQ_LEN_EST);
		if (seqs)
			pkg->block.seqs().reserve(record_est, data_.size());
		if (titles)
			pkg->block.reserve_ids(data_.size(), (BlockId)record_est);
		const char* ptr = data_.data(), *end = ptr + data_.size();
		int looks_like_dna = 0;
		while (ptr < end) {
			const char* record_end = find_record_end(ptr, end);
			if (*ptr != '>')
				throw std::runtime_error("Malformed FASTA record.");
			const char* title_end = std::find(ptr, record_end, '\n');
			if (titles)
				id.assign(ptr + 1, Util::String::trim_cr(ptr + 1, title_end));
			if (seqs)
				parse_sequence(title_end == record_end ? record_end : title_end + 1, record_end, seq);
			pkg->block.push_back(Sequence(seq), titles ? id.c_str() : nullptr, nullptr, std::numeric_limits<OId>::max(), seq_type, frame_mask);
			++pkg->seq_count;
			ptr = record_end == end ? end : record_end + 1;
			if (pkg->no == 0 && pkg->seq_count <= CHECK_FOR_DNA_COUNT && value_traits_.seq_type == SequenceType::amino_acid && Util::Seq::looks_like_dna(Sequence(seq)) && !config.ignore_warnings) {
				++looks_like_dna;
				if (looks_like_dna >= CHECK_FOR_DNA_COUNT)
					throw runtime_error(DNA_ERR);
			}
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
		if (end - begin <= 2)
			return end;
		const char* p = begin + 2;
		while (p < end) {
			const char* gt = static_cast<const char*>(std::memchr(p, '>', static_cast<size_t>(end - p)));
			if (!gt)
				return end;
			if (gt[-1] == '\n')
				return gt - 1;
			p = gt + 1;
		}
		return end;
	}

	void parse_sequence(const char* begin, const char* end, vector<Letter>& out) const {
		out.clear();
		out.reserve(end - begin);
		for (const char* ptr = begin; ptr < end;) {
			const char* line_end = std::find(ptr, end, '\n');
			const char* trimmed_end = Util::String::trim_cr(ptr, line_end);
			for (const char* i = ptr; i < trimmed_end; ++i)
				out.push_back(value_traits_.from_char(*i));
			ptr = line_end == end ? end : line_end + 1;
		}
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
	chunk->bytes_ += file_ptr_->read_to_fasta_record_end(chunk->data_);
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

optional<uint64_t> FastaFile::sequence_count() const {
	return seqs_ > 0 ? optional<uint64_t>(seqs_) : nullopt;
}

optional<uint64_t> FastaFile::letters() const {
	return letters_ > 0 ? optional<uint64_t>(letters_) : nullopt;
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
	return file_.front().name();
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

double FastaFile::rel_file_ptr()
{
	return (double)file_.front().tell() / file_.front().size();
}

pair<uint64_t, uint64_t> FastaFile::init_read() {
	uint64_t seqs = 0, letters = 0;
	const bool count_only = !flag_any(flags_, Flags::ACC_TO_OID_MAPPING | Flags::OID_TO_ACC_MAPPING | Flags::NEED_LENGTH_LOOKUP);
	if (index_file_.empty() || !exists(index_file_)) {
		const uint64_t limit = std::min<uint64_t>(Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT)), SequenceFile::DEFAULT_LOAD_SIZE);
		uint64_t bytes = 0, ms = 0;
		OId oid = 0;
		const bool map_seqids = flag_any(flags_, Flags::ACC_TO_OID_MAPPING | Flags::OID_TO_ACC_MAPPING);
		const SequenceFile::Flags flags = flags_;
		flags_ |= Flags::SEQS;
		if (map_seqids)
			flags_ |= Flags::TITLES;
		for (;;) {
			TaskTimer timer;
			Block* b = load_seqs(limit);
			*message_stream << std::setprecision(2) << std::fixed << rel_file_ptr() * 100 << "%" << endl;
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
			if (map_seqids) {
				seq_length_.reserve(seqs);
				for (BlockId i = 0; i < n; ++i) {
					add_seqid_mapping(b->ids()[i], oid);
					if (flag_any(flags_, Flags::NEED_LENGTH_LOOKUP))
						seq_length_.push_back(b->seqs().length(i));
					++oid;
				}
			}
			else if (flag_any(flags_, Flags::NEED_LENGTH_LOOKUP)) {
				const size_t offset = seq_length_.size();
				seq_length_.resize(offset + n);
				const size_t thread_count = std::min<size_t>(config.threads_, n);
				vector<std::thread> threads;
				threads.reserve(thread_count);
				for (size_t t = 0; t < thread_count; ++t) {
					const BlockId begin = (BlockId)(n * t / thread_count);
					const BlockId end = (BlockId)(n * (t + 1) / thread_count);
					threads.emplace_back([&, begin, end, offset]() {
						for (BlockId i = begin; i < end; ++i)
							seq_length_[offset + i] = b->seqs().length(i);
					});
				}
				for (std::thread& thread : threads)
					thread.join();
				oid += n;
			}
			delete b;
		}
		std::ostringstream ss;
		ss << "Loaded " << bytes << " bytes from disk at " << ((double)bytes / MEGABYTES / ms * 1e6) << " MB/s" << std::endl;
		open_stats_ += ss.str();
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
