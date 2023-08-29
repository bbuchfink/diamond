/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

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

#include <limits>
#include <fstream>
#include "../basic/config.h"
#include "../util/seq_file_format.h"
#include "../util/log_stream.h"
#include "../masking/masking.h"
#include "../taxonomy.h"
#include "../taxon_list.h"
#include "../taxonomy_nodes.h"
#include "../util/algo/MurmurHash3.h"
#include "../util/io/record_reader.h"
#include "../util/parallel/multiprocessing.h"
#include "dmnd.h"
#include "../reference.h"
#include "../taxonomy.h"
#include "../util/system/system.h"
#include "../util/algo/external_sort.h"
#include "../../util/util.h"
#include "../fasta/fasta_file.h"
#include "../../util/sequence/sequence.h"

using std::tuple;
using std::string;
using std::list;
using std::endl;
using std::unique_ptr;
using std::pair;
using std::runtime_error;

const char* DatabaseFile::FILE_EXTENSION = ".dmnd";
const uint32_t ReferenceHeader::current_db_version_prot = 3;
const uint32_t ReferenceHeader::current_db_version_nucl = 4;

Serializer& operator<<(Serializer &s, const ReferenceHeader2 &h)
{
	s.unset(Serializer::VARINT);
	s << sizeof(ReferenceHeader2);
	s.write(h.hash, sizeof(h.hash));
	s << h.taxon_array_offset << h.taxon_array_size << h.taxon_nodes_offset << h.taxon_names_offset;
#ifdef EXTRA
	s << (int32_t)h.db_type;
#endif
	return s;
}

Deserializer& operator>>(Deserializer &d, ReferenceHeader2 &h)
{
#ifdef EXTRA
	int32_t db_type;
#endif
	d.read_record().read(h.hash, sizeof(h.hash))
		>> h.taxon_array_offset
		>> h.taxon_array_size
		>> h.taxon_nodes_offset
		>> h.taxon_names_offset
#ifdef EXTRA
		>> db_type
#endif
		>> Finish();
#ifdef EXTRA
	h.db_type = (SequenceType)db_type;
#endif
	return d;
}

InputFile& operator>>(InputFile& file, ReferenceHeader& h)
{
	file.varint = false;
	file >> h.magic_number >> h.build >> h.db_version >> h.sequences >> h.letters >> h.pos_array_offset;
	return file;
}

Serializer& operator<<(Serializer& file, const ReferenceHeader& h)
{
	file.unset(Serializer::VARINT);
	file << h.magic_number << h.build << h.db_version << h.sequences << h.letters << h.pos_array_offset;
	return file;
}

InputFile& operator>>(InputFile& file, SequenceFile::SeqInfo& r) {
	uint32_t p;
	file >> r.pos >> r.seq_len >> p;
	return file;
}

Serializer& operator<<(Serializer& file, const SequenceFile::SeqInfo& r) {
	file << r.pos << r.seq_len << (uint32_t)0;
	return file;
}

SequenceFile::SeqInfo DatabaseFile::read_seqinfo() {
	SeqInfo r;
	(*this) >> r;
	pos_array_offset += SeqInfo::SIZE;
	return r;
}

void DatabaseFile::putback_seqinfo() {
	pos_array_offset -= SeqInfo::SIZE;
}

void DatabaseFile::init(Flags flags)
{
	read_header(*this, ref_header);
	if (flag_any(flags, Flags::NO_COMPATIBILITY_CHECK))
		return;
	if (ref_header.build < min_build_required || ref_header.db_version < MIN_DB_VERSION)
		throw std::runtime_error("Database was built with an older version of Diamond and is incompatible.");
	if (ref_header.db_version > std::max(ReferenceHeader::current_db_version_prot , ReferenceHeader::current_db_version_nucl))
		throw std::runtime_error("Database was built with a newer version of Diamond and is incompatible.");
	if (ref_header.sequences == 0)
		throw std::runtime_error("Incomplete database file. Database building did not complete successfully.");
	*this >> header2;
	pos_array_offset = ref_header.pos_array_offset;
}

DatabaseFile::DatabaseFile(const string &input_file, Metadata metadata, Flags flags, const ValueTraits& value_traits):
	SequenceFile(SequenceFile::Type::DMND, Alphabet::STD, flags, FormatFlags::DICT_LENGTHS | FormatFlags::DICT_SEQIDS | FormatFlags::SEEKABLE | FormatFlags::LENGTH_LOOKUP, value_traits),
	InputFile(auto_append_extension_if_exists(input_file, FILE_EXTENSION), InputFile::BUFFERED),
	temporary(false)
{
	init(flags);

	vector<string> e;
	if (flag_any(metadata, Metadata::TAXON_MAPPING) && !has_taxon_id_lists())
		e.push_back("taxonomy mapping information (--taxonmap option)");
	if (flag_any(metadata, Metadata::TAXON_NODES) && !has_taxon_nodes())
		e.push_back("taxonomy nodes information (--taxonnodes option)");
	if (flag_any(metadata, Metadata::TAXON_SCIENTIFIC_NAMES) && !has_taxon_scientific_names())
		e.push_back("taxonomy names information (--taxonnames option)");
	if (flag_any(metadata, Metadata::TAXON_RANKS) && build_version() < 131)
		e.push_back("taxonomy ranks information (database needs to be built with diamond version >= 0.9.30");

	if (!e.empty())
		throw std::runtime_error("Options require taxonomy information included in the database. Please use the respective options to build this information into the database when running diamond makedb: " + join(", ", e));

	if (flag_any(metadata, Metadata::TAXON_MAPPING))
		taxon_list_.reset(new TaxonList(seek(header2.taxon_array_offset), ref_header.sequences, header2.taxon_array_size));

	if (flag_any(metadata, Metadata::TAXON_SCIENTIFIC_NAMES)) {
		seek(header2.taxon_names_offset);
		(*this) >> taxon_scientific_names_;
	}

	if (flag_any(metadata, Metadata::TAXON_NODES))
		taxon_nodes_.reset(new TaxonomyNodes(seek(header2.taxon_nodes_offset), ref_header.build));

	if (flag_any(flags, Flags::ACC_TO_OID_MAPPING | Flags::OID_TO_ACC_MAPPING))
		read_seqid_list();
}

DatabaseFile::DatabaseFile(TempFile &tmp_file, const ValueTraits& value_traits):
	SequenceFile(SequenceFile::Type::DMND, Alphabet::STD, Flags::NONE, FormatFlags::DICT_LENGTHS | FormatFlags::DICT_SEQIDS | FormatFlags::SEEKABLE | FormatFlags::LENGTH_LOOKUP, value_traits),
	InputFile(tmp_file, 0),
	temporary(true)
{
	init();
}

int64_t DatabaseFile::file_count() const {
	return 1;
}

void DatabaseFile::close() {
	if (temporary)
		InputFile::close_and_delete();
	else
		InputFile::close();
}

void DatabaseFile::read_header(InputFile &stream, ReferenceHeader &header)
{
	stream >> header;
	if (header.magic_number != ReferenceHeader().magic_number)
		throw Database_format_exception();
}

bool DatabaseFile::has_taxon_id_lists() const
{
	return header2.taxon_array_offset != 0;
}

bool DatabaseFile::has_taxon_nodes() const
{
	return header2.taxon_nodes_offset != 0;
}

bool DatabaseFile::has_taxon_scientific_names() const {
	return header2.taxon_names_offset != 0;
}

static void push_seq(const Sequence &seq, const char *id, size_t id_len, uint64_t &offset, vector<SequenceFile::SeqInfo> &pos_array, OutputFile &out, size_t &letters, size_t &n_seqs)
{
	pos_array.emplace_back(offset, seq.length());
	out.write("\xff", 1);
	out.write(seq.data(), seq.length());
	out.write("\xff", 1);
	out.write(id, id_len + 1);
	letters += seq.length();
	++n_seqs;
	offset += seq.length() + id_len + 3;
}

void DatabaseFile::make_db()
{
	config.file_buffer_size = 4 * MEGABYTES;
	if (config.input_ref_file.size() > 1)
		throw std::runtime_error("Too many arguments provided for option --in.");
	const string input_file_name = config.input_ref_file.empty() ? string() : config.input_ref_file.front();
	if (input_file_name.empty())
		message_stream << "Input file parameter (--in) is missing. Input will be read from stdin." << endl;
	else
		message_stream << "Database input file: " << input_file_name << endl;

	TaskTimer total;
	TaskTimer timer("Opening the database file", true);

	value_traits = (config.dbtype == SequenceType::amino_acid) ? amino_acid_traits : nucleotide_traits;
    FastaFile db_file({ input_file_name }, Metadata (), Flags::NONE, value_traits);

    unique_ptr<OutputFile> out(new OutputFile(config.database));
	ReferenceHeader header;
    ReferenceHeader2 header2;

    *out << header;
	*out << header2;

	size_t letters = 0, n = 0, n_seqs = 0, total_seqs = 0;
	uint64_t offset = out->tell();

    LoadFlags flags = SequenceFile::LoadFlags::ALL;
    if (config.dbtype == SequenceType::nucleotide){
        header.db_version = ReferenceHeader::current_db_version_nucl;
        flags |= SequenceFile::LoadFlags::DNA_PRESERVATION;
    }

    Block* block;
	const FASTA_format format;
	vector<SeqInfo> pos_array;
	ExternalSorter<pair<string, OId>> accessions;
	AccessionParsing acc_stats;
	try {
		while (true) {
			timer.go("Loading sequences");
			block = db_file.load_seqs((int64_t)1e9, nullptr, flags);
			if (block->empty()) {
				delete block;
				break;
			}
			n = block->seqs().size();

			if (config.dbtype == SequenceType::amino_acid && config.masking_ != "0") {
				timer.go("Masking sequences");
				mask_seqs(block->seqs(), Masking::get(), false, MaskingAlgo::SEG);
			}

            timer.go("Writing sequences");
			for (size_t i = 0; i < n; ++i) {
				Sequence seq = block->seqs()[i];
				if (seq.length() == 0)
					throw std::runtime_error("File format error: sequence of length 0 at line " + std::to_string(db_file.line_count()));
				push_seq(seq, block->ids()[i], block->ids().length(i), offset, pos_array, *out, letters, n_seqs);
			}
			if (!config.prot_accession2taxid.empty()) {
				timer.go("Writing accessions");
				for (size_t i = 0; i < n; ++i) {
					vector<string> acc = accession_from_title(block->ids()[i], acc_stats);
					for (const string& s : acc)
						accessions.push(std::make_pair(s, total_seqs + i));
				}
			}
			timer.go("Hashing sequences");
			for (size_t i = 0; i < n; ++i) {
				Sequence seq = block->seqs()[i];
				MurmurHash3_x64_128(seq.data(), (int)seq.length(), header2.hash, header2.hash);
				MurmurHash3_x64_128(block->ids()[i], block->ids().length(i), header2.hash, header2.hash);
			}
			delete block;
			total_seqs += n;
		}
	}
	catch (std::exception&) {
		out->close();
		out->remove();
		throw;
	}

	timer.finish();

	timer.go("Writing trailer");
	header.pos_array_offset = offset;
	pos_array.emplace_back(offset, 0);
	for (const SeqInfo& r : pos_array)
		*out << r;
	pos_array.clear();
	pos_array.shrink_to_fit();
	timer.finish();

	if (!config.prot_accession2taxid.empty() && !config.no_parse_seqids)
		message_stream << endl << "Accession parsing rules triggered for database seqids (use --no-parse-seqids to disable):" << endl << acc_stats << endl;

	Util::Table stats;
	stats("Database sequences", n_seqs);
	stats("Database letters", letters);

	taxonomy.init();
	if (!config.prot_accession2taxid.empty()) {
		header2.taxon_array_offset = out->tell();
		TaxonList::build(*out, accessions, n_seqs, stats);
		header2.taxon_array_size = out->tell() - header2.taxon_array_offset;
	}
	if (!config.nodesdmp.empty()) {
		TaxonomyNodes nodes(config.nodesdmp);
		header2.taxon_nodes_offset = out->tell();
		nodes.save(*out);
	}
	if (!config.namesdmp.empty()) {
		header2.taxon_names_offset = out->tell();
		*out << taxonomy.name_;
	}

#ifdef EXTRA
    header2.db_type = config.dbtype;
#endif

    timer.go("Closing the input file");
	db_file.close();

	timer.go("Closing the database file");
	header.letters = letters;
	header.sequences = n_seqs;
	out->seek(0);
	*out << header;
	*out << header2;
	out->close();

	timer.finish();
	stats("Database hash", hex_print(header2.hash, 16));
	stats("Total time", total.get(), "s");

	message_stream << endl << stats;
}

void DatabaseFile::set_seqinfo_ptr(OId i) {
	pos_array_offset = ref_header.pos_array_offset + SeqInfo::SIZE * i;
}

OId DatabaseFile::tell_seq() const {
	return (pos_array_offset - ref_header.pos_array_offset) / SeqInfo::SIZE;
}

bool DatabaseFile::eof() const {
	return tell_seq() == sequence_count();
}

void DatabaseFile::init_seq_access() {
	seek(sizeof(ReferenceHeader) + sizeof(ReferenceHeader2) + 8);
}

bool DatabaseFile::read_seq(vector<Letter>& seq, string &id, std::vector<char>* quals)
{
	char c;
	read(&c, 1);
	seq.clear();
	id.clear();
	read_to(std::back_inserter(seq), '\xff');
	read_to(std::back_inserter(id), '\0');
	return false;
}

void DatabaseFile::skip_seq()
{
	char c;
	if(read(&c, 1) != 1)
		throw std::runtime_error("Unexpected end of file.");
	if(!seek_forward('\xff'))
		throw std::runtime_error("Unexpected end of file.");
	if(!seek_forward('\0'))
		throw std::runtime_error("Unexpected end of file.");
}

bool DatabaseFile::is_diamond_db(const string &file_name) {
	if (file_name == "-" || file_name.empty())
		return false;
	InputFile db_file(file_name);
	uint64_t magic_number = 0;
	try {
		db_file >> magic_number;
	}
	catch (EndOfStream&) {}
	bool r = (magic_number == ReferenceHeader::MAGIC_NUMBER);
	db_file.close();
	return r;
}

void DatabaseFile::create_partition_fixednumber(size_t n) {
	size_t max_letters_balanced = static_cast<size_t>(std::ceil(static_cast<double>(ref_header.letters)/static_cast<double>(n)));
	message_stream << "Fixed number partitioning using " << max_letters_balanced << " (" << n << ")" << endl;
	this->create_partition(max_letters_balanced);
}

void DatabaseFile::create_partition_balanced(size_t max_letters) {
	//double n = std::ceil(static_cast<double>(ref_header.letters) / static_cast<double>(max_letters));
	//size_t max_letters_balanced = static_cast<size_t>(std::ceil(static_cast<double>(ref_header.letters)/n));
	//cout << "Balanced partitioning using " << max_letters_balanced << " (" << max_letters << ")" << endl;
	this->create_partition(max_letters);
}

void DatabaseFile::create_partition(size_t max_letters) {
	TaskTimer timer("Create partition of DatabaseFile");
	size_t letters = 0, seqs = 0, total_seqs = 0;
	int i_chunk = 0;

	size_t oid = 0, oid_begin;
	set_seqinfo_ptr(oid);
	init_seqinfo_access();

	SeqInfo r = read_seqinfo();
	bool first = true;

	while (r.seq_len) {
		SeqInfo r_next = read_seqinfo();
		if (first) {
			oid_begin = oid;
			first = false;
		}
		letters += r.seq_len;
		++seqs;
		++total_seqs;
		if ((letters > max_letters) || (r_next.seq_len == 0)) {
			partition.chunks.push_back(Chunk(i_chunk, oid_begin, seqs));
			first = true;
			seqs = 0;
			letters = 0;
			++i_chunk;
		}
		r = r_next;
		++oid;
	}

	reverse(partition.chunks.begin(), partition.chunks.end());

	partition.max_letters = max_letters;
	partition.n_seqs_total = total_seqs;
}

int DatabaseFile::get_n_partition_chunks() {
	return (int)partition.chunks.size();
}

void DatabaseFile::save_partition(const string & partition_file_name, const string & annotation) {
	std::ofstream out(partition_file_name);
	// cout << "WRITING " << partition_file_name << endl;
	for (const auto& i : partition.chunks) {
		out << to_string(i);
		if (annotation.size() > 0) {
			out << " " << annotation;
		}
		out << endl;
	}
}

Chunk to_chunk(const string & line) {
	vector<string> tokens = split(line, ' ');
	return Chunk((int)stoull(tokens[0]), stoull(tokens[1]), stoull(tokens[2]));
}

string to_string(const Chunk & c) {
	const string buf = std::to_string(c.i) + " " + std::to_string(c.offset) + " " + std::to_string(c.n_seqs);
	return buf;
}

void DatabaseFile::load_partition(const string & partition_file_name) {
	string line;
	std::ifstream in(partition_file_name);
	clear_partition();
	while (getline(in, line)) {
		auto chunk = to_chunk(line);
		partition.chunks.push_back(chunk);
	}
}

void DatabaseFile::clear_partition() {
	partition.max_letters = 0;
	partition.n_seqs_total = 0;
	partition.chunks.clear();
}

void DatabaseFile::init_seqinfo_access() {
	seek(pos_array_offset);
}

void DatabaseFile::seek_chunk(const Chunk& chunk) {
	set_seqinfo_ptr(chunk.offset);
}

size_t DatabaseFile::id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) {
	return seq_info_next.pos - seq_info.pos - seq_info.seq_len - 3;
}

void DatabaseFile::seek_offset(size_t p) {
	seek(p);
}

void DatabaseFile::read_seq_data(Letter* dst, size_t len, size_t& pos, bool seek) {
	if (seek)
		this->seek(pos);
	read(dst - 1, len + 2);
	*(dst - 1) = Sequence::DELIMITER;
	*(dst + len) = Sequence::DELIMITER;
}

void DatabaseFile::read_id_data(const int64_t oid, char* dst, size_t len) {
	read(dst, len + 1);
}

void DatabaseFile::skip_id_data() {
	if (!seek_forward('\0')) throw std::runtime_error("Unexpected end of file.");
}

int64_t DatabaseFile::sequence_count() const {
	return ref_header.sequences;
}

size_t DatabaseFile::letters() const {
	return ref_header.letters;
}

int DatabaseFile::db_version() const {
	return ref_header.db_version;
}

int DatabaseFile::program_build_version() const {
	return ref_header.build;
}

SequenceFile::Metadata DatabaseFile::metadata() const {
	Metadata flags = Metadata();
	if (has_taxon_id_lists())
		flags |= Metadata::TAXON_MAPPING;
	if (has_taxon_nodes())
		flags |= Metadata::TAXON_NODES;
	if (has_taxon_scientific_names())
		flags |= Metadata::TAXON_SCIENTIFIC_NAMES;
	return flags;
}

int DatabaseFile::build_version() {
	return ref_header.build;
}

DatabaseFile::~DatabaseFile()
{
	close();
}

void DatabaseFile::close_weakly()
{
}

void DatabaseFile::reopen()
{
}

BitVector* DatabaseFile::filter_by_accession(const std::string& file_name)
{
	throw std::runtime_error("The .dmnd database format does not support filtering by accession.");
	return nullptr;
}

const BitVector* DatabaseFile::builtin_filter()
{
	return nullptr;
}

std::string DatabaseFile::file_name()
{
	return InputFile::file_name;
}

int64_t DatabaseFile::sparse_sequence_count() const
{
	return sequence_count();
}

std::vector<TaxId> DatabaseFile::taxids(size_t oid) const
{
	return (*taxon_list_)[oid];
}

void DatabaseFile::seq_data(size_t oid, std::vector<Letter>& dst) const
{
	throw OperationNotSupported();
}

size_t DatabaseFile::seq_length(size_t oid) const
{
	throw OperationNotSupported();
}

void DatabaseFile::init_random_access(const size_t query_block, const size_t ref_blocks, bool dictionary)
{
	if(dictionary)
		load_dictionary(query_block, ref_blocks);
}

void DatabaseFile::end_random_access(bool dictionary)
{
	if (!dictionary)
		return;
	free_dictionary();
}

void DatabaseFile::init_write() {
	throw OperationNotSupported();
}

void DatabaseFile::write_seq(const Sequence& seq, const std::string& id) {
	throw OperationNotSupported();
}

std::string DatabaseFile::taxon_scientific_name(TaxId taxid) const {
	return taxid < (TaxId)taxon_scientific_names_.size() && !taxon_scientific_names_[taxid].empty() ? taxon_scientific_names_[taxid] : std::to_string(taxid);
}

#ifdef EXTRA

const char* const SEQID_EXTENSION = ".seqid";

void DatabaseFile::prep_db() {
	config.database.require();
	TaskTimer timer("Opening the database");
	DatabaseFile db(config.database);
	timer.go("Writing seqid list");
	db.make_seqid_list(); // db.file_name() + SEQID_EXTENSION);
	timer.go("Closing the database");
	db.close();
}

#endif

void DatabaseFile::read_seqid_list() {
	if (flag_any(flags_, Flags::ACC_TO_OID_MAPPING))
		acc2oid_.reserve(sequence_count());
	OId oid = 0;
	vector<Letter> seq;
	string id;
	init_seq_access();
	for (int64_t n = 0; n < sequence_count(); ++n) {
		read_seq(seq, id);
		const char* msg = Util::Seq::fix_title(id);
		if (msg)
			message_stream << "Warning: " << msg << std::endl;
		add_seqid_mapping(id, oid++);
	}
	/*if ((flag_any(flags_, Flags::ACC_TO_OID_MAPPING) && (int64_t)acc2oid_.size() != sequence_count())
		|| (flag_any(flags_, Flags::OID_TO_ACC_MAPPING) && (int64_t)acc_.size() != sequence_count()))
		throw runtime_error("Inconsistent size of database and seqid file.");*/
}