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

#include <sqlite3.h>
#include <objmgr/object_manager.hpp>
#include <objmgr/scope.hpp>
#include <objmgr/util/create_defline.hpp>
#include <objtools/blast/blastdb_format/blastdb_dataextract.hpp>
#include <objtools/blast/seqdb_reader/impl/seqdbtax.hpp>
#include <corelib/ncbiutil.hpp>
#include <objtools/blast/seqdb_reader/seqdb.hpp>
//#include <objects/taxon1/local_taxon.hpp>
#include "util/io/text_input_file.h"
#include "blastdb.h"
#include "util/system/system.h"
#include "basic/config.h"
#include "util/log_stream.h"
#include "taxdmp.h"
#include "data/taxonomy_nodes.h"
//#include "phr.h"

using std::endl;
using std::vector;
using std::runtime_error;
using std::numeric_limits;
using namespace ncbi;

static string full_id(CBioseq& bioseq, CBioseq_Handle* bioseq_handle, bool long_ids, bool ctrl_a) {
	string id;
	if (long_ids) {
		CConstRef<CSeq_id> best_id = FindBestChoice(bioseq.GetId(), CSeq_id::FastaAARank);
		id = best_id->AsFastaString();
		sequence::CDeflineGenerator gen;
		id += gen.GenerateDefline(*bioseq_handle, 0);
	}
	else {
		CBlastDeflineUtil::ProcessFastaDeflines(bioseq, id, ctrl_a);
		id.erase(0, 1);
		id.pop_back();
	}
	return id;
}

template<typename It>
static void load_seq_data(CBioseq& bioseq, CBioseq_Handle bioseq_handle, It it) {
	ncbi::objects::CSeqVector v = bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac);
	if (v.GetCoding() != CSeq_data::e_Iupacaa)
		throw runtime_error("Invalid sequence coding in BLAST database.");
	
	for (ncbi::TSeqPos i = 0; i < v.size(); ++i) {
		const auto l = v[i] & 31;
		const Letter s = IUPACAA_TO_STD[l];
		if (s == -1)
			throw runtime_error("Unrecognized sequence character in BLAST database letter=" + std::to_string(l)
				+ " accession=" + bioseq.GetFirstId()->AsFastaString()
				+ " position=" + std::to_string(i + 1));
		*it = s;
		++it;
	}
}

list<CRef<CSeq_id>>::const_iterator best_id(const list<CRef<CSeq_id>>& ids) {
	if (ids.empty())
		throw runtime_error("Unable to retrieve sequence id from BLAST database.");
	auto min = ids.cbegin(), it = min;
	int min_score = (*min)->TextScore(), s;
	++it;
	while (it != ids.cend()) {
		if ((s = (*it)->TextScore()) < min_score) {
			min_score = s;
			min = it;
		}
		++it;
	}
	return min;
}

BlastDB::BlastDB(const std::string& file_name, Metadata metadata, Flags flags, const ValueTraits& value_traits) :
	SequenceFile(Type::BLAST, Alphabet::NCBI, flags, FormatFlags::TITLES_LAZY | FormatFlags::SEEKABLE | FormatFlags::LENGTH_LOOKUP | FormatFlags::DICT_LENGTHS | FormatFlags::DICT_SEQIDS, metadata, value_traits),
	file_name_(file_name),
	db_(new CSeqDB(file_name, CSeqDB::eProtein, nullptr, true)),
	taxon_db_(nullptr),
	oid_(0),
	long_seqids_(false),
	flags_(flags),
	sequence_count_(db_->GetNumOIDs()),
	sparse_sequence_count_(db_->GetNumSeqs())
{
	const string dbdir = containing_directory_absolute(SeqDB_ResolveDbPathNoExtension(db_->GetDBNameList(), 'p'));
	if (flag_any(metadata, Metadata::TAXON_RANKS)) {
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
	if (flag_any(metadata, Metadata::TAXON_NODES)) {
		const string dbpath = dbdir + PATH_SEPARATOR + "taxonomy4blast.sqlite3";
		if(!exists(dbpath))
			throw runtime_error("Taxonomy database file not found: " + dbpath + ". Make sure that the database was downloaded correctly.");
		/*CArgDescriptions desc;
		CLocalTaxon::AddArguments(desc);
		const char* argv[] = {"", "-taxon-db", dbpath.c_str(), nullptr};
		int argc = 3;
		CNcbiArguments ncbi_argv(argc, argv);
		unique_ptr<CArgs> args(desc.CreateArgs(ncbi_argv));
		taxon_.reset(new CLocalTaxon(*args));*/
		if (sqlite3_open_v2(dbpath.c_str(), &taxon_db_, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK) {
			string msg = taxon_db_ ? sqlite3_errmsg(taxon_db_) : "unknown error";
			if (taxon_db_) sqlite3_close(taxon_db_);
			throw runtime_error("Failed to open database " + dbpath + ": " + msg);
		}
		const TaxId max_id = max_taxid();
		parent_cache_.resize(max_id + 1, numeric_limits<TaxId>::min());
		init_cache();
	}
	if (flag_any(metadata, Metadata::TAXON_SCIENTIFIC_NAMES)) {
		if (getenv("BLASTDB") == nullptr)
#ifdef _MSC_VER
			_putenv_s("BLASTDB", dbdir.c_str());
#else
			if(setenv("BLASTDB", dbdir.c_str(), 0) != 0)
				cerr << "Warning: Failed to set BLASTDB environment variable: " << strerror(errno) << endl;
#endif
		SSeqDBTaxInfo info;
		if (!CSeqDBTaxInfo::GetTaxNames(1234, info))
			throw runtime_error("Taxonomy information (scientific names) is missing. Make sure that the BLAST database was downloaded correctly and that taxonomy files (taxdb.bti) are present in the database search paths (current working directory, BLAST database directory, BLASTDB environment variable)");
		const string file = dbdir + PATH_SEPARATOR + "names.dmp";
		if (!exists(file))
			throw runtime_error("Taxonomy names information (names.dmp) is missing in search path ("
				+ dbdir + "). Download and extract this file in the database directory: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip");
		auto f = [&](TaxId taxid, const string& name) {			
			if (!CSeqDBTaxInfo::GetTaxNames(taxid, info))
				extra_names_[taxid] = name;
			};
		read_names_dmp(file, f);
	}
	if (config.multiprocessing)
		throw runtime_error("Multiprocessing mode is not compatible with BLAST databases.");
}

void BlastDB::print_info() const {
	if (flag_any(metadata_, Metadata::TAXON_RANKS)) {
		if (!custom_ranks_.empty())
			message_stream << "Custom taxonomic ranks in database: " << custom_ranks_.size() << endl;
		message_stream << "Taxonomic ids assigned to ranks: " << rank_mapping_.size() << endl;
	}
	if (flag_any(metadata_, Metadata::TAXON_NODES))
		message_stream << "Maximum taxid in database: " << parent_cache_.size() - 1 << endl;
	if (flag_any(metadata_, Metadata::TAXON_SCIENTIFIC_NAMES))
		message_stream << "Extra taxonomic scientific names in names.dmp: " << extra_names_.size() << endl;
}

int64_t BlastDB::file_count() const {
	return 1;
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
	if (oid_ >= db_->GetNumOIDs()) {
		++oid_;
		return SeqInfo(0, 0);
	}
	const int l = db_->GetSeqLength(oid_);
	if (l == 0)
		throw std::runtime_error("Database with sequence length 0 is not supported");
	return SeqInfo(oid_++, l);
}

void BlastDB::putback_seqinfo()
{
	--oid_;
}

size_t BlastDB::id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next)
{
	if(flag_any(flags_, Flags::FULL_TITLES))
		return full_id(*db_->GetBioseq((BlastOid)seq_info.pos), nullptr, long_seqids_, true).length();
	else {
		return (*best_id(db_->GetSeqIDs((BlastOid)seq_info.pos)))->GetSeqIdString(true).length();
	}
}

void BlastDB::seek_offset(size_t p)
{
}

void BlastDB::read_seq_data(Letter* dst, size_t len, size_t& pos, bool seek)
{
	*(dst - 1) = Sequence::DELIMITER;
	*(dst + len) = Sequence::DELIMITER;
	const char* buf;
	const int db_len = db_->GetSequence((int)pos, &buf);
	std::copy(buf, buf + len, dst);
	++pos;
	db_->RetSequence(&buf);
}

void BlastDB::read_id_data(const int64_t oid, char* dst, size_t len)
{
	if (flag_any(flags_, Flags::FULL_TITLES)) {
		const string id = full_id(*db_->GetBioseq((BlastOid)oid), nullptr, long_seqids_, true);
		std::copy(id.begin(), id.begin() + len, dst);
	}
	else {
		const string id = (*best_id(db_->GetSeqIDs((BlastOid)oid)))->GetSeqIdString(true);
		std::copy(id.begin(), id.end(), dst);
	}
	dst[len] = '\0';
}

void BlastDB::skip_id_data()
{
	//++oid_seqdata_;
}

string BlastDB::fetch_seqid(OId oid, bool all) const {
	if (flag_any(flags_, Flags::FULL_TITLES)) {
		return full_id(*db_->GetBioseqNoData((int)oid), nullptr, long_seqids_, true);
	}
	else {
		const list<CRef<CSeq_id>> ids = db_->GetSeqIDs(oid);
		if (ids.empty())
			return "N/A";
		if (all) {
			string r;
			for (list<CRef<CSeq_id>>::const_iterator i = ids.begin(); i != ids.end(); ++i) {
				if (i != ids.begin())
					r += '\1';
				r += (*i)->GetSeqIdString(true);
			}
			return r;
		}
		else
			return ids.front()->GetSeqIdString(true);
	}
}

string BlastDB::seqid(OId oid, bool all) const
{
	return fetch_seqid(oid, all);
}

Loc BlastDB::dict_len(DictId dict_id, const size_t ref_block) const
{
	if (dict_id >= (int64_t)dict_oid_[dict_block(ref_block)].size())
		throw std::runtime_error("Dictionary not loaded.");
	return db_->GetSeqLength((int)dict_oid_[dict_block(ref_block)][dict_id]);
}

vector<Letter> BlastDB::dict_seq(DictId dict_id, const size_t ref_block) const
{
	const size_t b = dict_block(ref_block);
	if (dict_id >= (DictId)dict_oid_[b].size())
		throw std::runtime_error("Dictionary not loaded.");
	vector<Letter> v;
	seq_data(dict_oid_[b][dict_id], v);
	alph_ncbi_to_std(v.begin(), v.end());
	return v;
}

int64_t BlastDB::sequence_count() const
{
	return sequence_count_;
}

int64_t BlastDB::sparse_sequence_count() const
{
	return sparse_sequence_count_;
}

int64_t BlastDB::letters() const
{
	return db_->GetTotalLength();
}

int BlastDB::db_version() const
{
	return (int)db_->GetBlastDbVersion();
}

int BlastDB::program_build_version() const
{
	return 0;
}

bool BlastDB::read_seq(vector<Letter>& seq, string& id, vector<char>* quals)
{
	id.clear();
	CRef<CBioseq> bioseq = db_->GetBioseq(oid_);
	CScope scope(*CObjectManager::GetInstance());
	CBioseq_Handle bioseq_handle = scope.AddBioseq(*bioseq);

	id = full_id(*bioseq, &bioseq_handle, long_seqids_, false);

	seq.clear();
	load_seq_data(*bioseq, bioseq_handle, std::back_inserter(seq));
	
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

void BlastDB::set_seqinfo_ptr(int64_t i)
{
	oid_ = (int)i;
}

void BlastDB::close()
{
	db_.reset();
	sqlite3_close(taxon_db_);
}

BitVector* BlastDB::filter_by_accession(const string& file_name)
{
	BitVector* v = new BitVector(sequence_count());
	TextInputFile in(file_name);
	vector<string> accs;
	while (in.getline(), (!in.line.empty() || !in.eof())) {
		accs.push_back(in.line);
	}
	in.close();

	vector<CSeqDB::TOID> oids;
	try {
		db_->AccessionsToOids(accs, oids);
	}
	catch (CSeqDBException& e) {
		throw std::runtime_error(e.GetMsg());
	}

	for (size_t i = 0; i < accs.size(); ++i) {
		if (oids[i] < 0)
			if (config.skip_missing_seqids)
				message_stream << "WARNING: Accession not found in database : " + accs[i] << endl;
			else
				throw std::runtime_error("Accession not found in database: " + accs[i] + ". Use --skip-missing-seqids to ignore.");
		else
			v->set(oids[i]);
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
	db_->GetTaxIDs((BlastOid)oid, v);
	return v;
}

TaxId BlastDB::max_taxid() const {
	set<TaxId> taxids;
	db_->GetDBTaxIds(taxids);
	TaxId m = taxids.empty() ? 0 : *taxids.rbegin();
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
			m = std::max(m, sqlite3_column_int(stmt, 0));
		} else {
			string msg = sqlite3_errmsg(taxon_db_);
			sqlite3_finalize(stmt);
			sqlite3_close(taxon_db_);
			throw runtime_error("SQLite step error: " + msg);
		}
	}
	return m;
}

TaxId BlastDB::get_parent(TaxId taxid) {
	//return taxon_->GetParent(taxid);
	if (taxid <= 0)
		return taxid;
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
	SSeqDBTaxInfo info;
	if(CSeqDBTaxInfo::GetTaxNames(taxid, info))
		return info.scientific_name;
	else {
		auto it = extra_names_.find(taxid);
		return it == extra_names_.end() ? std::to_string(taxid) : it->second;
	}
}

void BlastDB::seq_data(size_t oid, vector<Letter>& dst) const
{
	const char* buf;
	const int db_len = db_->GetSequence((int)oid, &buf);
	dst.clear();
	dst.resize(db_len);
	std::copy(buf, buf + db_len, dst.data());
	db_->RetSequence(&buf);
}

size_t BlastDB::seq_length(size_t oid) const
{
	return db_->GetSeqLength((int)oid);
}

void BlastDB::end_random_access(bool dictionary)
{
	if(dictionary)
		free_dictionary();
}

vector<OId> BlastDB::accession_to_oid(const string& acc) const
{
	vector<int> r;
	db_->AccessionToOids(acc, r);
	if(r.empty())
		throw runtime_error("Accession not found in database: " + acc);
	return vector<OId>(r.begin(), r.end());
}

const BitVector* BlastDB::builtin_filter() {
	if (sequence_count() == sparse_sequence_count())
		return nullptr;
	if (oid_filter_.empty()) {
		oid_filter_ = BitVector(sequence_count());
		int oid = 0;
		while (db_->CheckOrFindOID(oid)) {
			oid_filter_.set(oid);
			++oid;
		}
	}
	return &oid_filter_;
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

StringSet BlastDB::load_ids(OId begin, OId end) const {
	return StringSet();
	/*auto it = volumes_.lower_bound(begin);
	if(it == volumes_.end())
		throw runtime_error("No volume found for OID " + std::to_string(begin));
	Phr phr(it->second + ".phr", it->second + ".pin");
	StringSet ids;
	for (OId i = 0; i < phr.size(); ++i) {
		ids.reserve(phr.len(i));
	}
	ids.finish_reserve();
	vector<pair<const uint8_t*, size_t>> titles;
	vector<string> seqids;
	vector<uint64_t> taxids;
	for (OId i = 0; i < phr.size(); ++i) {
		if(!phr.parse_record(i, titles, seqids, taxids))
			throw runtime_error("Failed to parse PHR record with OID " + std::to_string(i) + " in BLAST database.");
		ids.assign(i, seqids[0].begin(), seqids[0].end());
	}
	return ids;*/
}