/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
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
#include "../basic/config.h"
#include "../util/seq_file_format.h"
#include "../util/log_stream.h"
#include "../basic/masking.h"
#include "../taxonomy.h"
#include "../util/io/file_backed_buffer.h"
#include "../taxon_list.h"
#include "../taxonomy_nodes.h"
#include "../util/algo/MurmurHash3.h"
#include "../util/io/record_reader.h"
#include "../util/parallel/multiprocessing.h"
#include "dmnd.h"
#include "../reference.h"
#include "../load_seqs.h"
#include "../taxonomy.h"

using namespace std;

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

void DatabaseFile::init(int flags)
{
	read_header(*this, ref_header);
	if (flags & NO_COMPATIBILITY_CHECK)
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

DatabaseFile::DatabaseFile(const string &input_file, int flags):
	SequenceFile(SequenceFile::Type::DMND),
	InputFile(input_file, InputFile::BUFFERED),
	temporary(false)
{
	init(flags);
}

DatabaseFile::DatabaseFile(TempFile &tmp_file):
	SequenceFile(SequenceFile::Type::DMND),
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

	size_t letters = 0, n = 0, n_seqs = 0;
	uint64_t offset = out->tell();

	SequenceSet *seqs;
	String_set<char, 0> *ids;
	const FASTA_format format;
	vector<SeqInfo> pos_array;
	FileBackedBuffer accessions;

	try {
		while ((timer.go("Loading sequences"), n = ::load_seqs(db_file->begin(), db_file->end(), format, &seqs, ids, 0, nullptr, (size_t)(1e9), string(), amino_acid_traits)) > 0) {
			if (config.masking == 1) {
				timer.go("Masking sequences");
				mask_seqs(*seqs, Masking::get(), false);
			}
			timer.go("Writing sequences");
			for (size_t i = 0; i < n; ++i) {
				Sequence seq = (*seqs)[i];
				if (seq.length() == 0)
					throw std::runtime_error("File format error: sequence of length 0 at line " + to_string(db_file->front().line_count));
				push_seq(seq, (*ids)[i], ids->length(i), offset, pos_array, *out, letters, n_seqs);
			}
			if (!config.prot_accession2taxid.empty()) {
				timer.go("Writing accessions");
				for (size_t i = 0; i < n; ++i)
					accessions << Taxonomy::Accession::from_title((*ids)[i]);
			}
			timer.go("Hashing sequences");
			for (size_t i = 0; i < n; ++i) {
				Sequence seq = (*seqs)[i];
				MurmurHash3_x64_128(seq.data(), (int)seq.length(), header2.hash, header2.hash);
				MurmurHash3_x64_128((*ids)[i], ids->length(i), header2.hash, header2.hash);
			}
			delete seqs;
			delete ids;
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
	timer.finish();

	taxonomy.init();
	if (!config.prot_accession2taxid.empty()) {
		header2.taxon_array_offset = out->tell();
		TaxonList::build(*out, accessions.rewind(), n_seqs);
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
	message_stream << "Database hash = " << hex_print(header2.hash, 16) << endl;
	message_stream << "Processed " << n_seqs << " sequences, " << letters << " letters." << endl;
	message_stream << "Total time = " << total.get() << "s" << endl;
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
	double n = std::ceil(static_cast<double>(ref_header.letters) / static_cast<double>(max_letters));
	size_t max_letters_balanced = static_cast<size_t>(std::ceil(static_cast<double>(ref_header.letters)/n));
	cout << "Balanced partitioning using " << max_letters_balanced << " (" << max_letters << ")" << endl;
	this->create_partition(max_letters_balanced);
}

void DatabaseFile::create_partition(size_t max_letters) {
	task_timer timer("Create partition of DatabaseFile");
	size_t letters = 0, seqs = 0, total_seqs = 0;
	size_t i_chunk = 0;

	rewind();
	seek(pos_array_offset);

	SeqInfo r, r_next;
	size_t pos;
	bool first = true;

	read(&r, 1);
	while (r.seq_len) {
		if (first) {
			pos = pos_array_offset;
			first = false;
		}
		letters += r.seq_len;
		++seqs;
		++total_seqs;
		read(&r_next, 1);
		if ((letters > max_letters) || (r_next.seq_len == 0)) {
			partition.chunks.push_back(Chunk(i_chunk, pos, seqs));
			first = true;
			seqs = 0;
			letters = 0;
			++i_chunk;
		}
		pos_array_offset += SeqInfo::SIZE;
		r = r_next;
	}

	reverse(partition.chunks.begin(), partition.chunks.end());

	partition.max_letters = max_letters;
	partition.n_seqs_total = total_seqs;
}

size_t DatabaseFile::get_n_partition_chunks() {
	return partition.chunks.size();
}

void DatabaseFile::save_partition(const string & partition_file_name, const string & annotation) {
	ofstream out(partition_file_name);
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
	const string buf = to_string(c.i) + " " + to_string(c.offset) + " " + to_string(c.n_seqs);
	return buf;
}

void DatabaseFile::load_partition(const string & partition_file_name) {
	string line;
	ifstream in(partition_file_name);
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
	seek(chunk.offset);
}

size_t DatabaseFile::id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next) {
	return seq_info_next.pos - seq_info.pos - seq_info.seq_len - 3;
}

void DatabaseFile::seek_offset(size_t p) {
	seek(p);
}

void DatabaseFile::read_seq_data(Letter* dst, size_t len) {
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

void DatabaseFile::check_metadata(int flags) const {
	if ((flags & TAXON_MAPPING) && !has_taxon_id_lists())
		throw runtime_error("Output format requires taxonomy mapping information built into the database (use --taxonmap parameter for the makedb command).");
	if ((flags & TAXON_NODES) && !has_taxon_nodes())
		throw runtime_error("Output format requires taxonomy nodes information built into the database (use --taxonnodes parameter for the makedb command).");
	if ((flags & TAXON_SCIENTIFIC_NAMES) && !has_taxon_scientific_names())
		throw runtime_error("Output format requires taxonomy names information built into the database (use --taxonnames parameter for the makedb command).");
}

int DatabaseFile::metadata() const {
	int flags = 0;
	if (has_taxon_id_lists())
		flags |= TAXON_MAPPING;
	if (has_taxon_nodes())
		flags |= TAXON_NODES;
	if (has_taxon_scientific_names())
		flags |= TAXON_SCIENTIFIC_NAMES;
	return flags;
}

TaxonList* DatabaseFile::taxon_list() {
	return new TaxonList(seek(header2.taxon_array_offset), ref_header.sequences, header2.taxon_array_size);
}

TaxonomyNodes* DatabaseFile::taxon_nodes() {
	return new TaxonomyNodes(seek(header2.taxon_nodes_offset), ref_header.build);
}

int DatabaseFile::build_version() {
	return ref_header.build;
}

std::vector<string>* DatabaseFile::taxon_scientific_names() {
	vector<string>* r = new vector<string>;
	seek(header2.taxon_names_offset);
	(*this) >> (*r);
	return r;
}