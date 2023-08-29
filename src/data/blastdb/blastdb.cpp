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

#include <objmgr/object_manager.hpp>
#include <objmgr/scope.hpp>
#include <objmgr/util/create_defline.hpp>
#include <objtools/blast/blastdb_format/blastdb_dataextract.hpp>
#include <objtools/blast/seqdb_reader/impl/seqdbtax.hpp>
#include <corelib/ncbiutil.hpp>
#include "../util/io/text_input_file.h"
#include "blastdb.h"
#include "../../util/string/tokenizer.h"
#include "../../util/system/system.h"
#include "../basic/config.h"
#include "../../util/util.h"

using std::cout;
using std::endl;
using std::vector;
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

template<typename _it>
static void load_seq_data(CBioseq& bioseq, CBioseq_Handle bioseq_handle, _it it) {
	ncbi::objects::CSeqVector v = bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac);
	if (v.GetCoding() != CSeq_data::e_Iupacaa)
		throw std::runtime_error("Invalid sequence coding in BLAST database.");
	
	for (ncbi::TSeqPos i = 0; i < v.size(); ++i) {
		const auto l = v[i] & 31;
		const Letter s = IUPACAA_TO_STD[l];
		if (s == -1)
			throw std::runtime_error("Unrecognized sequence character in BLAST database letter=" + std::to_string(l)
				+ " accession=" + bioseq.GetFirstId()->AsFastaString()
				+ " position=" + std::to_string(i + 1));
		*it = s;
		++it;
	}
}

list<CRef<CSeq_id>>::const_iterator best_id(const list<CRef<CSeq_id>>& ids) {
	if (ids.empty())
		throw std::runtime_error("Unable to retrieve sequence id from BLAST database.");
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
	SequenceFile(Type::BLAST, Alphabet::NCBI, flags, FormatFlags::TITLES_LAZY | FormatFlags::SEEKABLE | FormatFlags::LENGTH_LOOKUP, value_traits),
	file_name_(file_name),
	db_(new CSeqDBExpert(file_name, CSeqDB::eProtein)),
	oid_(0),
	long_seqids_(false),
	flags_(flags),
	sequence_count_(db_->GetNumOIDs()),
	sparse_sequence_count_(db_->GetNumSeqs())
{
#ifndef EXTRA
	if (flag_any(metadata, Metadata::TAXON_NODES | Metadata::TAXON_MAPPING | Metadata::TAXON_SCIENTIFIC_NAMES | Metadata::TAXON_RANKS))
		throw std::runtime_error("Taxonomy features are not supported for the BLAST database format.");
#endif
	vector<string> paths;
	CSeqDB::FindVolumePaths(file_name, CSeqDB::eProtein, paths);
	for (const string& db : paths)
		if(!exists(db + ".acc"))
			throw std::runtime_error("Accession file not found. BLAST databases require preprocessing using this command line: diamond prepdb -d DATABASE_PATH");
	if (config.multiprocessing)
		throw std::runtime_error("Multiprocessing mode is not compatible with BLAST databases.");

	if (flag_any(metadata, Metadata::TAXON_NODES)) {
		const string path = extract_dir(SeqDB_ResolveDbPathNoExtension(db_->GetDBNameList(), 'p'));
		if (path.empty())
			throw std::runtime_error("Could not find BLAST database path.");
		try {
			taxon_nodes_.reset(new TaxonomyNodes(path + dir_separator + "nodes.dmp", true));
		}
		catch (FileOpenException&) {
			throw std::runtime_error("Taxonomy nodes file (nodes.dmp) was not found in the BLAST database directory.");
		}
	}
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

std::string BlastDB::seqid(OId oid) const
{
	if (flag_any(flags_, Flags::FULL_TITLES)) {
		return full_id(*db_->GetBioseq((int)oid), nullptr, long_seqids_, true);
	}
	else {
		if (oid >= acc_.size())
			throw std::runtime_error("Accession array not correctly initialized.");
		return acc_[oid];
	}
}

std::string BlastDB::dict_title(DictId dict_id, const size_t ref_block) const
{
	if (dict_id >= (DictId)dict_oid_[dict_block(ref_block)].size())
		throw std::runtime_error("Dictionary not loaded.");
	return seqid(dict_oid_[dict_block(ref_block)][dict_id]);
}

Loc BlastDB::dict_len(DictId dict_id, const size_t ref_block) const
{
	if (dict_id >= (int64_t)dict_oid_[dict_block(ref_block)].size())
		throw std::runtime_error("Dictionary not loaded.");
	return db_->GetSeqLength((int)dict_oid_[dict_block(ref_block)][dict_id]);
}

std::vector<Letter> BlastDB::dict_seq(DictId dict_id, const size_t ref_block) const
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

size_t BlastDB::letters() const
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

bool BlastDB::read_seq(std::vector<Letter>& seq, std::string& id, std::vector<char>* quals)
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

SequenceFile::Metadata BlastDB::metadata() const
{
	return Metadata();
}

int BlastDB::build_version()
{
	return 0;
}

void BlastDB::create_partition_balanced(size_t max_letters)
{
	throw OperationNotSupported();
}

void BlastDB::save_partition(const std::string& partition_file_name, const std::string& annotation)
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
}

void BlastDB::close_weakly()
{
	db_.reset();
}

void BlastDB::reopen()
{
	if(db_.get() == nullptr)
		db_.reset(new CSeqDBExpert(file_name_, CSeqDB::eProtein));
}

BitVector* BlastDB::filter_by_accession(const std::string& file_name)
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

std::string BlastDB::file_name()
{
	return file_name_;
}

std::vector<TaxId> BlastDB::taxids(size_t oid) const
{
	vector<TaxId> v;
	db_->GetTaxIDs((BlastOid)oid, v);
	return v;
}

string BlastDB::taxon_scientific_name(TaxId taxid) const {
	SSeqDBTaxInfo info;
	if(CSeqDBTaxInfo::GetTaxNames(taxid, info))
		return info.scientific_name;
	else
		return std::to_string(taxid);
}

void BlastDB::seq_data(size_t oid, std::vector<Letter>& dst) const
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

const char* BlastDB::ACCESSION_FIELD = "#accession*";

void BlastDB::init_random_access(const size_t query_block, const size_t ref_block, bool dictionary)
{
	reopen();
	if(dictionary)
		load_dictionary(query_block, ref_block);
	if (flag_any(flags_, Flags::FULL_TITLES))
		return;
	TaskTimer timer("Loading accessions");
	vector<string> paths;
	if(db_.get())
		db_->FindVolumePaths(paths);
	else
		CSeqDB::FindVolumePaths(file_name_, CSeqDB::eProtein, paths);
	acc_.clear();
	string acc;
	for (const string& path : paths) {
		TextInputFile f(path + ".acc");
		f.getline();
		if (f.line != ACCESSION_FIELD)
			throw std::runtime_error("Accession file is missing the header field: " + path);
		while (f.getline(), !f.eof() || !f.line.empty()) {
			if (flag_any(flags_, Flags::ALL_SEQIDS))
				acc = Util::String::replace(f.line, '\t', '\1');
			else
				Util::String::Tokenizer(f.line, "\t") >> acc;
			acc_.push_back(acc.begin(), acc.end());
		}
		f.close();
	}
}

void BlastDB::end_random_access(bool dictionary)
{
	acc_.clear();
	acc_.shrink_to_fit();
	if(dictionary)
		free_dictionary();
}

std::vector<OId> BlastDB::accession_to_oid(const std::string& acc) const
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

void BlastDB::prep_blast_db(const string& path) {
	vector<string> paths;
	CSeqDB::FindVolumePaths(path, CSeqDB::eProtein, paths);
	for (const string& db : paths) {
		message_stream << "Processing volume: " << db << endl;
		CSeqDB volume(db, CSeqDB::eProtein);
		const int n = volume.GetNumOIDs();
		message_stream << "Number of sequences: " << n << endl;
		OutputFile out(db + ".acc", Compressor::ZSTD);
		TextBuffer buf;
		buf << BlastDB::ACCESSION_FIELD << '\n';
		size_t id_count = 0;
		for (int i = 0; i < n; ++i) {
			list<CRef<CSeq_id>> ids = volume.GetSeqIDs(i);
			if (!ids.empty()) {
				//auto best = best_id(ids);
				auto it = ids.cbegin();
				buf << (*it)->GetSeqIdString(true);
				while (++it != ids.cend()) {
					//if (it != best)
						buf << '\t' << (*it)->GetSeqIdString(true);
				}
				id_count += ids.size();
			}
			buf << '\n';
			out.write(buf.data(), buf.size());
			buf.clear();
		}
		message_stream << "Number of accessions: " << id_count << endl;
		out.close();
	}
}

void BlastDB::init_write() {
	throw OperationNotSupported();
}

void BlastDB::write_seq(const Sequence& seq, const std::string& id) {
	throw OperationNotSupported();
}
