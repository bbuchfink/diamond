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
#include "../basic/masking.h"
#include "../taxonomy.h"
#include "../taxon_list.h"
#include "../taxonomy_nodes.h"
#include "../util/algo/MurmurHash3.h"
#include "../util/io/record_reader.h"
#include "../util/parallel/multiprocessing.h"
#include "dmnd.h"
#include "../reference.h"
#include "../load_seqs.h"
#include "../taxonomy.h"
#include "../util/system/system.h"
#include "../util/algo/external_sort.h"
#include "../../util/util.h"

using std::tuple;
using std::string;
using std::list;
using std::cout;
using std::endl;

const char* DatabaseFile::FILE_EXTENSION = ".dmnd";

Serializer& operator<<(Serializer &s, const ReferenceHeader2 &h)
{
	s.unset(Serializer::VARINT);
	s << sizeof(ReferenceHeader2);
	s.write(h.hash, sizeof(h.hash));
	s << h.taxon_array_offset << h.taxon_array_size << h.taxon_nodes_offset << h.taxon_names_offset;
	return s;
}

Deserializer& operator>>(Deserializer &d, ReferenceHeader2 &h)
{
	d.read_record().read(h.hash, sizeof(h.hash))
		>> h.taxon_array_offset
		>> h.taxon_array_size
		>> h.taxon_nodes_offset
		>> h.taxon_names_offset
		>> Finish();
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

InputFile& operator>>(InputFile& file, SeqInfo& r) {
	uint32_t p;
	file >> r.pos >> r.seq_len >> p;
	return file;
}

Serializer& operator<<(Serializer& file, const SeqInfo& r) {
	file << r.pos << r.seq_len << (uint32_t)0;
	return file;
}

SeqInfo DatabaseFile::read_seqinfo() {
	SeqInfo r;
	(*this) >> r;
	pos_array_offset += SeqInfo::SIZE;
	return r;
}

void DatabaseFile::putback_seqinfo() {
	pos_array_offset -= SeqInfo::SIZE;
}

void DatabaseFile::write_dict_entry(size_t block, size_t oid, size_t len, const char* id, const Letter* seq)
{
	OutputFile& f = *dict_file_;
	f << (uint32_t)oid;
	f << (uint32_t)len;
	f << id;
	if (flag_any(flags_, Flags::TARGET_SEQS))
		f.write(seq, len);
	dict_alloc_size_ += strlen(id);
}

bool DatabaseFile::load_dict_entry(InputFile& f, size_t ref_block)
{
	uint32_t oid, len;
	string title;
	try {
		f >> oid;
	}
	catch (EndOfStream&) {
		return false;
	}
	f >> len >> title;
	const size_t b = dict_block(ref_block);
	dict_oid_[b].push_back(oid);
	dict_len_[b].push_back(len);
	dict_title_[b].push_back(title.begin(), title.end());
	if (flag_any(flags_, Flags::TARGET_SEQS)) {
		vector<Letter> v(len);
		f.read(v.data(), len);
		dict_seq_[b].push_back(v.begin(), v.end());
	}
	return true;
}

void DatabaseFile::reserve_dict(const size_t ref_blocks)
{
	if (config.multiprocessing) {
		dict_len_ = std::vector<std::vector<uint32_t>>(ref_blocks);
		dict_title_ = std::vector<StringSet>(ref_blocks);
		if (flag_any(flags_, Flags::TARGET_SEQS))
			dict_seq_ = std::vector<SequenceSet>(ref_blocks);
	}
	else {
		dict_len_ = { {} };
		dict_len_[0].reserve(next_dict_id_);
		dict_title_ = { {} };
		dict_title_[0].reserve(next_dict_id_, dict_alloc_size_);
		if (flag_any(flags_, Flags::TARGET_SEQS)) {
			dict_seq_ = { {} };
			dict_seq_[0].reserve(next_dict_id_, 0);
		}
	}
}

void DatabaseFile::init(Flags flags)
{
	read_header(*this, ref_header);
	if (flag_any(flags, Flags::NO_COMPATIBILITY_CHECK))
		return;
	if (ref_header.build < min_build_required || ref_header.db_version < MIN_DB_VERSION)
		throw std::runtime_error("Database was built with an older version of Diamond and is incompatible.");
	if (ref_header.db_version > ReferenceHeader::current_db_version)
		throw std::runtime_error("Database was built with a newer version of Diamond and is incompatible.");
	if (ref_header.sequences == 0)
		throw std::runtime_error("Incomplete database file. Database building did not complete successfully.");
	*this >> header2;
	pos_array_offset = ref_header.pos_array_offset;
}

DatabaseFile::DatabaseFile(const string &input_file, Metadata metadata, Flags flags):
	SequenceFile(SequenceFile::Type::DMND, Alphabet::STD, flags),
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
}

DatabaseFile::DatabaseFile(TempFile &tmp_file):
	SequenceFile(SequenceFile::Type::DMND, Alphabet::STD, Flags::NONE),
	InputFile(tmp_file, 0),
	temporary(true)
{
	init();
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

static void push_seq(const Sequence &seq, const char *id, size_t id_len, uint64_t &offset, vector<SeqInfo> &pos_array, OutputFile &out, size_t &letters, size_t &n_seqs)
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

void DatabaseFile::make_db(TempFile **tmp_out, list<TextInputFile> *input_file)
{
	config.file_buffer_size = 4 * MEGABYTES;
	if (config.input_ref_file.size() > 1)
		throw std::runtime_error("Too many arguments provided for option --in.");
	const string input_file_name = config.input_ref_file.empty() ? string() : config.input_ref_file.front();
	if (input_file_name.empty() && !input_file)
		std::cerr << "Input file parameter (--in) is missing. Input will be read from stdin." << endl;
	if(!input_file && !input_file_name.empty())
		message_stream << "Database input file: " << input_file_name << endl;

	task_timer total;
	task_timer timer("Opening the database file", true);
	list<TextInputFile>* db_file;
	if (input_file)
		db_file = input_file;
	else {
		db_file = new list<TextInputFile>;
		db_file->emplace_back(input_file_name);
	}

	OutputFile *out = tmp_out ? new TempFile() : new OutputFile(config.database);
	ReferenceHeader header;
	ReferenceHeader2 header2;

	*out << header;
	*out << header2;

	size_t letters = 0, n = 0, n_seqs = 0, total_seqs = 0;
	uint64_t offset = out->tell();

	Block* block;
	const FASTA_format format;
	vector<SeqInfo> pos_array;
	ExternalSorter<pair<string, uint32_t>> accessions;

	try {
		while (true) {
			timer.go("Loading sequences");
			block = new Block(db_file->begin(), db_file->end(), format, (size_t)(1e9), amino_acid_traits, false);
			if (block->empty()) {
				delete block;
				break;
			}
			n = block->seqs().size();
			if (config.masking == 1) {
				timer.go("Masking sequences");
				mask_seqs(block->seqs(), Masking::get(), false);
			}
			timer.go("Writing sequences");
			for (size_t i = 0; i < n; ++i) {
				Sequence seq = block->seqs()[i];
				if (seq.length() == 0)
					throw std::runtime_error("File format error: sequence of length 0 at line " + std::to_string(db_file->front().line_count));
				push_seq(seq, block->ids()[i], block->ids().length(i), offset, pos_array, *out, letters, n_seqs);
			}
			if (!config.prot_accession2taxid.empty()) {
				timer.go("Writing accessions");
				for (size_t i = 0; i < n; ++i) {
					vector<string> acc = accession_from_title(block->ids()[i]);
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

	Table stats;
	stats("Database sequences", n_seqs);
	stats("Database letters", letters);
	taxonomy.init();
	if (!config.prot_accession2taxid.empty()) {
		header2.taxon_array_offset = out->tell();
		TaxonList::build(*out, accessions, n_seqs, stats);
		header2.taxon_array_size = out->tell() - header2.taxon_array_offset;
	}
	if (!config.nodesdmp.empty()) {
		header2.taxon_nodes_offset = out->tell();
		TaxonomyNodes::build(*out);
	}
	if (!config.namesdmp.empty()) {
		header2.taxon_names_offset = out->tell();
		*out << taxonomy.name_;
	}

	if (!input_file) {
		timer.go("Closing the input file");
		db_file->front().close();
		delete db_file;
	}

	timer.go("Closing the database file");
	header.letters = letters;
	header.sequences = n_seqs;
	out->seek(0);
	*out << header;
	*out << header2;
	if (tmp_out) {
		*tmp_out = static_cast<TempFile*>(out);
	} else {
		out->close();
		delete out;
	}

	timer.finish();
	stats("Database hash", hex_print(header2.hash, 16));
	stats("Total time", total.get(), "s");

	message_stream << endl << stats;
}

void DatabaseFile::set_seqinfo_ptr(size_t i) {
	pos_array_offset = ref_header.pos_array_offset + SeqInfo::SIZE * i;
}

size_t DatabaseFile::tell_seq() const {
	return (pos_array_offset - ref_header.pos_array_offset) / SeqInfo::SIZE;
}

void DatabaseFile::init_seq_access() {
	seek(sizeof(ReferenceHeader) + sizeof(ReferenceHeader2) + 8);
}

void DatabaseFile::read_seq(vector<Letter>& seq, string &id)
{
	char c;
	read(&c, 1);
	seq.clear();
	id.clear();
	read_to(std::back_inserter(seq), '\xff');
	read_to(std::back_inserter(id), '\0');
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
	if (file_name == "-")
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
	cout << "Fixed number partitioning using " << max_letters_balanced << " (" << n << ")" << endl;
	this->create_partition(max_letters_balanced);
}

void DatabaseFile::create_partition_balanced(size_t max_letters) {
	//double n = std::ceil(static_cast<double>(ref_header.letters) / static_cast<double>(max_letters));
	//size_t max_letters_balanced = static_cast<size_t>(std::ceil(static_cast<double>(ref_header.letters)/n));
	//cout << "Balanced partitioning using " << max_letters_balanced << " (" << max_letters << ")" << endl;
	this->create_partition(max_letters);
}

void DatabaseFile::create_partition(size_t max_letters) {
	task_timer timer("Create partition of DatabaseFile");
	size_t letters = 0, seqs = 0, total_seqs = 0;
	size_t i_chunk = 0;

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

size_t DatabaseFile::get_n_partition_chunks() {
	return partition.chunks.size();
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
	return Chunk(stoull(tokens[0]), stoull(tokens[1]), stoull(tokens[2]));
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
	current_ref_block = chunk.i;
	set_seqinfo_ptr(chunk.offset);
}

std::string DatabaseFile::seqid(size_t oid) const
{
	throw std::runtime_error("Operation not supported (seqid).");
}
std::string DatabaseFile::dict_title(size_t dict_id, const size_t ref_block) const
{
	const size_t b = dict_block(ref_block);
	if (b >= dict_title_.size() || dict_id >= dict_title_[b].size())
		throw std::runtime_error("Dictionary not loaded.");
	return dict_title_[b][dict_id];
}
size_t DatabaseFile::dict_len(size_t dict_id, const size_t ref_block) const
{
	const size_t b = dict_block(ref_block);
	if (b >= dict_len_.size() || dict_id >= dict_len_[b].size())
		throw std::runtime_error("Dictionary not loaded.");
	return dict_len_[b][dict_id];
}

std::vector<Letter> DatabaseFile::dict_seq(size_t dict_id, const size_t ref_block) const
{
	const size_t b = dict_block(ref_block);
	if (b >= dict_seq_.size() || dict_id >= dict_seq_[b].size())
		throw std::runtime_error("Dictionary not loaded.");
	Sequence s = dict_seq_[b][dict_id];
	return vector<Letter>(s.data(), s.end());
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

void DatabaseFile::read_id_data(char* dst, size_t len) {
	read(dst, len + 1);
}

void DatabaseFile::skip_id_data() {
	if (!seek_forward('\0')) throw std::runtime_error("Unexpected end of file.");
}

size_t DatabaseFile::sequence_count() const {
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

TaxonomyNodes* DatabaseFile::taxon_nodes() {
	return new TaxonomyNodes(seek(header2.taxon_nodes_offset), ref_header.build);
}

int DatabaseFile::build_version() {
	return ref_header.build;
}

DatabaseFile::~DatabaseFile()
{
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

BitVector* DatabaseFile::filter_by_taxonomy(const std::string& include, const std::string& exclude, TaxonomyNodes& nodes)
{
	if (!taxon_list_.get())
		throw std::runtime_error("Database does not contain taxonomy mapping.");
	BitVector* v = new BitVector(taxon_list_->size());
	if (!include.empty() && !exclude.empty())
		throw std::runtime_error("Options --taxonlist and --taxon-exclude are mutually exclusive.");
	const bool e = !exclude.empty();
	const std::set<unsigned> taxon_filter_list(parse_csv(e ? exclude : include));
	if (taxon_filter_list.empty())
		throw std::runtime_error("Option --taxonlist/--taxon-exclude used with empty list.");
	if (taxon_filter_list.find(1) != taxon_filter_list.end() || taxon_filter_list.find(0) != taxon_filter_list.end())
		throw std::runtime_error("Option --taxonlist/--taxon-exclude used with invalid argument (0 or 1).");
	for (size_t i = 0; i < taxon_list_->size(); ++i)
		if (nodes.contained((*taxon_list_)[i], taxon_filter_list) ^ e)
			v->set(i);
	return v;
}

const BitVector* DatabaseFile::builtin_filter()
{
	return nullptr;
}

std::string DatabaseFile::file_name()
{
	return InputFile::file_name;
}

size_t DatabaseFile::sparse_sequence_count() const
{
	return sequence_count();
}

std::vector<unsigned> DatabaseFile::taxids(size_t oid) const
{
	return (*taxon_list_)[oid];
}

void DatabaseFile::seq_data(size_t oid, std::vector<Letter>& dst) const
{
	throw std::runtime_error("Operation not supported.");
}

size_t DatabaseFile::seq_length(size_t oid) const
{
	throw std::runtime_error("Operation not supported.");
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
	dict_len_.clear();
	dict_len_.shrink_to_fit();
	dict_title_.clear();
	dict_title_.shrink_to_fit();
	dict_seq_.clear();
	dict_seq_.shrink_to_fit();
}

SequenceFile::LoadTitles DatabaseFile::load_titles()
{
	return LoadTitles::SINGLE_PASS;
}

std::vector<string>* DatabaseFile::taxon_scientific_names() {
	vector<string>* r = new vector<string>;
	seek(header2.taxon_names_offset);
	(*this) >> (*r);
	return r;
}