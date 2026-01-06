/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <sqlite3.h>
#include <iostream>
#include "util/io/text_input_file.h"
#include "blastdb.h"
#include "util/system/system.h"
#include "basic/config.h"
#include "util/log_stream.h"
#include "taxdmp.h"
#include "data/taxonomy_nodes.h"
#include "volume.h"

using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::runtime_error;
using std::numeric_limits;
using std::set;
using std::pair;

BlastDB::BlastDB(const string& file_name, Flags flags, const ValueTraits& value_traits) :
	SequenceFile(Type::BLAST, flags, FormatFlags::SEEKABLE | FormatFlags::LENGTH_LOOKUP | FormatFlags::DICT_LENGTHS | FormatFlags::DICT_SEQIDS, value_traits),
	file_name_(file_name),
	pal_(file_name),
	taxon_db_(nullptr),
	oid_(0),
	long_seqids_(false),
	flags_(flags),
	volume_(pal_.volumes[0], 0, pal_.oid_index[0], pal_.oid_index[1], true),
	raw_chunk_no_(0)
{
	if (pal_.metadata.find("SEQIDLIST") != pal_.metadata.end()) {
		flags_ |= Flags::NEED_LENGTH_LOOKUP;
	}
	if (pal_.metadata.find("TAXIDLIST") != pal_.metadata.end()) {
		flags_ |= SequenceFile::Flags::TAXON_MAPPING;
		flags_ |= SequenceFile::Flags::TAXON_NODES;
		flags_ |= SequenceFile::Flags::NEED_EARLY_TAXON_MAPPING | SequenceFile::Flags::NEED_LENGTH_LOOKUP;
	}
	if (bool(flags_ & Flags::TAXON_MAPPING)) {
		taxon_mapping_.reserve(pal_.sequence_count);
		//if (bool(flags_ & Flags::NEED_EARLY_TAXON_MAPPING)) {
		const auto flags_now = flags_;
		flags_ &= ~(Flags::SEQS | Flags::TITLES);
		for (;;) {
			pair<Block*, int64_t> b = load_parallel(1000000000, nullptr, nullptr, Chunk(), true);
			delete b.first;
			if (b.second == 0)
				break;
		}
		flags_ = flags_now;
		flags_ &= ~Flags::TAXON_MAPPING;
	}
	if (bool(flags_ & Flags::NEED_LENGTH_LOOKUP)) {
		seq_length_.reserve(pal_.sequence_count);
		for (OId i = 0; i < pal_.sequence_count; ++i) {
			open_volume(i);
			const uint32_t volume_oid = uint32_t(i - volume_.begin);
			seq_length_.push_back(volume_.length(volume_oid));
		}
	}
	string dbdir, dbfile;
	tie(dbdir, dbfile) = absolute_path(file_name);
	if (flag_any(flags_, Flags::TAXON_RANKS)) {
		const string file = dbdir + PATH_SEPARATOR + "nodes.dmp";
		if(!exists(file))
			throw runtime_error("Taxonomy rank information (nodes.dmp) is missing in search path ("
				+ dbdir + "). Download and extract this file in the database directory: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip");
		int next_rank = Rank::count;
		auto f = [&](TaxId taxid, TaxId parent, const string& rank) {
			const int p = Rank::predefined(rank.c_str());
			if (p >= 0) {
				rank_mapping_[taxid] = p;
			}
			else {
				auto it = custom_ranks_.find(rank);
				if (it == custom_ranks_.end()) {
					const int r = next_rank++;
					custom_ranks_[rank] = r;
					rank_mapping_[taxid] = r;
				}
				else
					rank_mapping_[taxid] = it->second;
			}
			};
		read_nodes_dmp(file, f);
	}
	if (flag_any(flags_, Flags::TAXON_NODES)) {
		const string dbpath = dbdir + PATH_SEPARATOR + "taxonomy4blast.sqlite3";
		if(!exists(dbpath))
			throw runtime_error("Taxonomy database (taxonomy4blast.sqlite3) file not found in path: " + dbpath + ". Make sure that the database was downloaded correctly.");
		if (sqlite3_open_v2(dbpath.c_str(), &taxon_db_, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK) {
			string msg = taxon_db_ ? sqlite3_errmsg(taxon_db_) : "unknown error";
			if (taxon_db_) sqlite3_close(taxon_db_);
			throw runtime_error("Failed to open database " + dbpath + ": " + msg);
		}
		const TaxId max_id = max_taxid();
		parent_cache_.resize(max_id + 1, numeric_limits<TaxId>::min());
		init_cache();
	}
	if (flag_any(flags_, Flags::TAXON_SCIENTIFIC_NAMES)) {
		//throw runtime_error("Taxonomy information (scientific names) is missing. Make sure that the BLAST database was downloaded correctly and that taxonomy files (taxdb.bti) are present in the database search paths (current working directory, BLAST database directory, BLASTDB environment variable)");
		const string file = dbdir + PATH_SEPARATOR + "names.dmp";
		if (!exists(file))
			throw runtime_error("Taxonomy names information (names.dmp) is missing in search path ("
				+ dbdir + "). Download and extract this file in the database directory: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip");
		auto f = [&](TaxId taxid, const string& name) {
			if (extra_names_.find(taxid) == extra_names_.end())
				extra_names_.emplace(taxid, name);
			};
		read_names_dmp(file, f);
	}
	if (config.multiprocessing)
		throw runtime_error("Multiprocessing mode is not compatible with BLAST databases.");
}

int64_t BlastDB::file_count() const {
	return 1;
}

void BlastDB::open_volume(OId oid) {
	if (oid >= volume_.begin && oid < volume_.end)
		return;
	const int idx = pal_.volume(oid);
	volume_ = Volume(pal_.volumes[idx], idx, pal_.oid_index[idx], pal_.oid_index[idx + 1], true);
}

void BlastDB::print_info() const {
	message_stream << "Database: " << config.database << ' ';
	message_stream << "(type: BLAST database, ";
	message_stream << "volumes: " << pal_.volumes.size() << ", ";
	message_stream << "sequences: " << sequence_count() << ", ";
	message_stream << "letters: " << letters() << ')' << endl;
	if (flag_any(flags_, Flags::TAXON_RANKS)) {
		if (!custom_ranks_.empty())
			message_stream << "Custom taxonomic ranks in database: " << custom_ranks_.size() << endl;
		message_stream << "Taxonomic ids assigned to ranks: " << rank_mapping_.size() << endl;
	}
	if (flag_any(flags_, Flags::TAXON_NODES))
		message_stream << "Maximum taxid in database: " << parent_cache_.size() - 1 << endl;
	if (flag_any(flags_, Flags::TAXON_SCIENTIFIC_NAMES))
		message_stream << "Extra taxonomic scientific names in names.dmp: " << extra_names_.size() << endl;
}

void BlastDB::init_seqinfo_access()
{
}

void BlastDB::init_seq_access()
{
	oid_ = 0;
}

void BlastDB::seek_chunk(const Chunk& chunk)
{
	throw OperationNotSupported();
}

OId BlastDB::tell_seq() const
{
	return oid_;
}

bool BlastDB::eof() const {
	return oid_ == sequence_count();
}

SequenceFile::SeqInfo BlastDB::read_seqinfo()
{
	if (oid_ >= pal_.sequence_count) {
		++oid_;
		return SeqInfo(0, 0);
	}
	const Loc l = seq_length(oid_);
	if (l == 0)
		throw runtime_error("Database with sequence length 0 is not supported");
	return SeqInfo(oid_++, l);
}

void BlastDB::putback_seqinfo()
{
	--oid_;
}

size_t BlastDB::id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next)
{
	open_volume(seq_info.pos);
	const uint32_t volume_oid = uint32_t(seq_info.pos - volume_.begin);
	return volume_.id_len(volume_oid);
}

void BlastDB::seek_offset(size_t p)
{
}

RawChunk* BlastDB::raw_chunk(size_t letters, SequenceFile::Flags flags) {
	Volume::RawChunk* c = volume_.raw_chunk(letters, flags);
	c->no = raw_chunk_no_++;
	OId oid = volume_.begin + volume_.seq_ptr();
	if (oid < pal_.sequence_count)
		open_volume(oid);
	return c;
}

void BlastDB::read_seq_data(Letter* dst, size_t len, size_t& pos, bool seek)
{
	*(dst - 1) = Sequence::DELIMITER;
	*(dst + len) = Sequence::DELIMITER;
	open_volume(pos);
	const uint32_t volume_oid = uint32_t(pos - volume_.begin);
	const vector<Letter> seq = volume_.sequence(volume_oid);
	copy(seq.begin(), seq.end(), dst);
	++pos;
}

void BlastDB::read_id_data(const int64_t oid, char* dst, size_t len, bool all, bool full_titles)
{
	const string id = fetch_seqid(oid, all, full_titles);
	copy(id.begin(), id.end(), dst);
	dst[len] = '\0';
}

vector<BlastDefLine> BlastDB::deflines(OId oid, bool all, bool full_titles, bool taxids) {
	open_volume(oid);
	const uint32_t volume_oid = uint32_t(oid - volume_.begin);
	return volume_.deflines(volume_oid, all, full_titles, taxids);
}

void BlastDB::skip_id_data()
{
	//++oid_seqdata_;
}

string BlastDB::fetch_seqid(OId oid, bool all, bool full_titles) {
	open_volume(oid);
	bool taxids = flag_any(flags_, Flags::TAXON_MAPPING);
	const uint32_t volume_oid = uint32_t(oid - volume_.begin);
	const vector<BlastDefLine> deflines = volume_.deflines(volume_oid, all, full_titles, taxids);
	if (taxids && taxon_mapping_.find(oid) == taxon_mapping_.end()) {
		for (auto i = deflines.begin(); i != deflines.end(); ++i)
			if (i->taxid)
				taxon_mapping_.emplace(oid, i->taxid.value());
	}
	return build_title(deflines, "\1", all);
}

void BlastDB::add_taxid_mapping(const std::vector<std::pair<OId, TaxId>>& taxids) {
	for(const auto& t : taxids) {
		taxon_mapping_.emplace(t.first, t.second);
	}
}

string BlastDB::seqid(OId oid, bool all, bool full_titles)
{
	return fetch_seqid(oid, all, full_titles);
}

/*Loc BlastDB::dict_len(DictId dict_id, const size_t ref_block)
{
	if (dict_id >= (int64_t)dict_oid_[dict_block(ref_block)].size())
		throw runtime_error("Dictionary not loaded.");
	return seq_length((int)dict_oid_[dict_block(ref_block)][dict_id]);
}*/

vector<Letter> BlastDB::dict_seq(DictId dict_id, const size_t ref_block)
{
	const size_t b = dict_block(ref_block);
	if (dict_id >= (DictId)dict_oid_[b].size())
		throw runtime_error("Dictionary not loaded.");
	vector<Letter> v;
	seq_data(dict_oid_[b][dict_id], v);
	return v;
}

uint64_t BlastDB::sequence_count() const
{
	return pal_.sequence_count;
}

uint64_t BlastDB::letters() const
{
	return pal_.letters;
}

int BlastDB::db_version() const
{
	return pal_.version;
}

int BlastDB::program_build_version() const
{
	return 0;
}

bool BlastDB::read_seq(vector<Letter>& seq, string& id, vector<char>* quals)
{
	open_volume(oid_);
	const uint32_t volume_oid = uint32_t(oid_ - volume_.begin);
	seq = volume_.sequence(volume_oid);
	const vector<BlastDefLine> deflines = volume_.deflines(volume_oid, true, true, false);
	id = build_title(deflines, " >", true);
	++oid_;
	return true;
}

int BlastDB::build_version()
{
	return 0;
}

void BlastDB::create_partition_balanced(int64_t max_letters)
{
	throw OperationNotSupported();
}

void BlastDB::save_partition(const string& partition_file_name, const string& annotation)
{
	throw OperationNotSupported();
}

int BlastDB::get_n_partition_chunks()
{
	throw OperationNotSupported();
}

void BlastDB::set_seqinfo_ptr(OId i)
{
	if(i != 0)
		throw runtime_error("Setting seqinfo pointer to non-zero value is not supported in BLAST databases.");
	oid_ = i;
	raw_chunk_no_ = 0;
	if (volume_.begin == 0)
		volume_.rewind();
	else
		open_volume(0);
}

void BlastDB::close()
{
	sqlite3_close(taxon_db_);
}

DbFilter* BlastDB::filter_by_accession(const string& file_name)
{
	DbFilter* v = new DbFilter(sequence_count());
	TextInputFile in(file_name);
	std::unordered_map<string, bool> accs;
	while (in.getline(), (!in.line.empty() || !in.eof())) {
		accs.emplace(in.line, false);
	}
	in.close();

	for (;;) {
		pair<Block*, int64_t> b = load_parallel(1000000000, nullptr, &accs, Chunk(), false);
		const BlockId n = b.first->oid_count();
		for (BlockId i = 0; i < n; ++i) {
			const OId oid = b.first->block_id2oid(i);
			v->oid_filter.set(oid);
			v->letter_count += seq_length(oid);
		}
		delete b.first;
		if (b.second == 0)
			break;
	}

	if (!config.skip_missing_seqids) {
		for (const auto& a : accs) {
			if(!a.second)
				throw runtime_error("Accession not found in database: " + a.first + ". Use --skip-missing-seqids to ignore.");
		}
	}
	return v;
}

string BlastDB::file_name()
{
	return file_name_;
}

vector<TaxId> BlastDB::taxids(size_t oid) const
{
	vector<TaxId> v;
	auto range = taxon_mapping_.equal_range(oid);
	for(auto it = range.first; it != range.second; ++it)
		v.push_back(it->second);
	return v;
}

TaxId BlastDB::max_taxid() const {
	if (taxon_db_) {
		static const char* sql = "SELECT max(taxid) FROM TaxidInfo;";
		sqlite3_stmt* stmt = nullptr;
		if (sqlite3_prepare_v2(taxon_db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
			string msg = sqlite3_errmsg(taxon_db_);
			sqlite3_close(taxon_db_);
			throw runtime_error("Failed to prepare statement: " + msg);
		}
		int rc = sqlite3_step(stmt);
		if (rc == SQLITE_ROW) {
			return sqlite3_column_int(stmt, 0);
		} else {
			string msg = sqlite3_errmsg(taxon_db_);
			sqlite3_finalize(stmt);
			sqlite3_close(taxon_db_);
			throw runtime_error("SQLite step error: " + msg);
		}
	}
	else
		return 0;
}

TaxId BlastDB::get_parent(TaxId taxid) {
	if (taxid <= 0)
		return taxid;
	if (taxid >= parent_cache_.size())
		return -1;
	if (parent_cache_[taxid] != numeric_limits<TaxId>::min())
		return parent_cache_[taxid];
	static const char* sql = "SELECT parent FROM TaxidInfo WHERE taxid = ?1 LIMIT 1;";

	sqlite3_stmt* stmt = nullptr;
	if (sqlite3_prepare_v2(taxon_db_, sql, -1, &stmt, nullptr) != SQLITE_OK) {
		string msg = sqlite3_errmsg(taxon_db_);
		sqlite3_close(taxon_db_);
		throw runtime_error("Failed to prepare statement: " + msg);
	}

	if (sqlite3_bind_int(stmt, 1, taxid) != SQLITE_OK) {
		string msg = sqlite3_errmsg(taxon_db_);
		sqlite3_finalize(stmt);
		sqlite3_close(taxon_db_);
		throw runtime_error("Failed to bind parameter: " + msg);
	}

	TaxId result;
	int rc = sqlite3_step(stmt);
	if (rc == SQLITE_ROW) {
		result = sqlite3_column_int(stmt, 0);
	}
	else if (rc == SQLITE_DONE) {
		result = -1;
	}
	else {
		string msg = sqlite3_errmsg(taxon_db_);
		sqlite3_finalize(stmt);
		sqlite3_close(taxon_db_);
		throw runtime_error("SQLite step error: " + msg);
	}

	parent_cache_[taxid] = result;
	sqlite3_finalize(stmt);
	return result;
}

string BlastDB::taxon_scientific_name(TaxId taxid) const {
	auto it = extra_names_.find(taxid);
	return it == extra_names_.end() ? std::to_string(taxid) : it->second;
}

void BlastDB::seq_data(size_t oid, vector<Letter>& dst)
{
	open_volume(oid);
	const uint32_t volume_oid = uint32_t(oid_ - volume_.begin);
	dst = volume_.sequence(volume_oid);
}

Loc BlastDB::seq_length(size_t oid) {
	if (oid < seq_length_.size()) {
		return seq_length_[oid];
	}
	else {
		open_volume(oid);
		const uint32_t volume_oid = uint32_t(oid - volume_.begin);
		return volume_.length(volume_oid);
	}
}

void BlastDB::end_random_access(bool dictionary)
{
	if(dictionary)
		free_dictionary();
}

vector<OId> BlastDB::accession_to_oid(const string& acc) const
{
	vector<int> r;
	//db_->AccessionToOids(acc, r);
	if(r.empty())
		throw runtime_error("Accession not found in database: " + acc);
	return vector<OId>(r.begin(), r.end());
}

BlastDB::~BlastDB()
{
	close();
}

void BlastDB::init_write() {
	throw OperationNotSupported();
}

void BlastDB::write_seq(const Sequence& seq, const string& id) {
	throw OperationNotSupported();
}