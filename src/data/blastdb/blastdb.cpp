#include <objmgr/object_manager.hpp>
#include <objmgr/scope.hpp>
#include <objmgr/util/create_defline.hpp>
#include <objtools/blast/blastdb_format/blastdb_dataextract.hpp>
#include <corelib/ncbiutil.hpp>
#include "blastdb.h"

using std::cout;
using std::endl;
using std::vector;
using namespace ncbi;

BlastDB::BlastDB(const std::string& file_name) :
	SequenceFile(Type::BLAST),
	db_(file_name, ncbi::CSeqDB::eProtein),
	oid_(0),
	long_seqids_(false)
{
}

void BlastDB::init_seqinfo_access()
{
}

void BlastDB::init_seq_access()
{
}

void BlastDB::seek_chunk(const Chunk& chunk)
{
}

size_t BlastDB::tell_seq() const
{
	return size_t();
}

SeqInfo BlastDB::read_seqinfo()
{
	return SeqInfo();
}

void BlastDB::putback_seqinfo()
{
}

size_t BlastDB::id_len(const SeqInfo& seq_info, const SeqInfo& seq_info_next)
{
	return size_t();
}

void BlastDB::seek_offset(size_t p)
{
}

void BlastDB::read_seq_data(Letter* dst, size_t len)
{
}

void BlastDB::read_id_data(char* dst, size_t len)
{
}

void BlastDB::skip_id_data()
{
}

size_t BlastDB::sequence_count() const
{
	return db_.GetNumSeqs();
}

size_t BlastDB::letters() const
{
	return db_.GetTotalLength();
}

int BlastDB::db_version() const
{
	return (int)db_.GetBlastDbVersion();
}

int BlastDB::program_build_version() const
{
	return 0;
}

void BlastDB::read_seq(std::vector<Letter>& seq, std::string& id)
{
	id.clear();
	CRef<CBioseq> bioseq = db_.GetBioseq(oid_);
	CScope scope(*CObjectManager::GetInstance());
	CBioseq_Handle bioseq_handle = scope.AddBioseq(*bioseq);

	if (long_seqids_) {
		CConstRef<CSeq_id> best_id = FindBestChoice(bioseq->GetId(), CSeq_id::FastaAARank);
		id = best_id->AsFastaString();
		sequence::CDeflineGenerator gen;
		id += gen.GenerateDefline(bioseq_handle, 0);
	}
	else {
		CBlastDeflineUtil::ProcessFastaDeflines(*bioseq, id, false);
		id.erase(0, 1);
		id.pop_back();
	}
	
	ncbi::objects::CSeqVector v = bioseq_handle.GetSeqVector(CBioseq_Handle::eCoding_Iupac);
	if (v.GetCoding() != CSeq_data::e_Iupacaa)
		throw std::runtime_error("Invalid sequence coding in BLAST database.");

	seq.clear();
	seq.reserve(v.size());
	for (size_t i = 0; i < v.size(); ++i) {
		const auto l = v[i] & 31;
		const Letter s = IUPACAA_TO_STD[l];
		if (s == -1)
			throw std::runtime_error("Unrecognized sequence character in BLAST database Letter=" + std::to_string(l)
				+ " Accession=" + bioseq->GetFirstId()->AsFastaString()
				+ " position=" + std::to_string(i + 1));
		seq.push_back(s);
	}
	++oid_;
}

void BlastDB::check_metadata(int flags) const
{
}

int BlastDB::metadata() const
{
	return 0;
}

TaxonList* BlastDB::taxon_list()
{
	return nullptr;
}

TaxonomyNodes* BlastDB::taxon_nodes()
{
	return nullptr;
}

std::vector<string>* BlastDB::taxon_scientific_names()
{
	return nullptr;
}

int BlastDB::build_version()
{
	return 0;
}

void BlastDB::create_partition_balanced(size_t max_letters)
{
}

void BlastDB::save_partition(const std::string& partition_file_name, const std::string& annotation)
{
}

size_t BlastDB::get_n_partition_chunks()
{
	return size_t();
}

void BlastDB::set_seqinfo_ptr(size_t i)
{
}

void BlastDB::close()
{
}
