#include <fstream>
#include "common.h"
#include "../../output/output_format.h"
#include "../cluster.h"

using std::endl;

namespace Cluster { namespace Incremental {

Config::Config() :
	message_stream(true),
	verbosity(1),
	sens(Search::cluster_sens.at(config.sensitivity)),
	block_size(config.chunk_size == 0.0 ? Search::sensitivity_traits[0].at(config.sensitivity).default_block_size : config.chunk_size),
	output_format(init_output(-1)),
	db(SequenceFile::auto_create({ config.database }, SequenceFile::Flags::NEED_LETTER_COUNT | SequenceFile::Flags::OID_TO_ACC_MAPPING)),
	centroids(config.resume.empty() ? new FastaFile(config.output_file + ".centroids.faa", true, FastaFile::WriteAccess())
	: new FastaFile(config.resume + ".centroids.faa", false, FastaFile::WriteAccess())),
	output_file(open_out_tsv()),
	seqs_processed(0),
	letters_processed(0),
	oid2centroid(db->sequence_count()),
	time_self_aln(0),
	time_search(sens.size() + 1, 0),
	problem_size(sens.size() + 1, 0),
	problem_size_self(0)
{
	sens.push_back(config.sensitivity);
	for (int i = 0; i < (int)sens.size() - 1; ++i)
		cache.emplace_back(new Block);
}

Config::~Config() {
}

void Config::status_msg() {
	message_stream << "#Seqs=" << seqs_processed << " #Centroids=" << centroids->sequence_count() << " Time=" << total_time.seconds() << "s" << std::endl;
}

void Config::save_state() {
	std::ofstream out1(config.output_file + ".oid2centroid");
	for (OId i = 0; i < db->tell_seq(); ++i)
		out1 << oid2centroid[i] << endl;
	std::ofstream out2(config.output_file + ".centroid2oid");
	for (OId i : centroid2oid)
		out2 << i << endl;
}

void Config::load_state() {
	std::ifstream in1(config.resume + ".oid2centroid");
	int64_t n = 0;
	CentroidId i;
	while (in1 >> i)
		oid2centroid[n++] = i;
	std::ifstream in2(config.resume + ".centroid2oid");
	OId oid;
	while (in2 >> oid)
		centroid2oid.push_back(oid);
	this->message_stream << "Centroid count = " << centroid2oid.size() << endl;
	this->message_stream << "Seeking to OId " << n << endl;
	db->set_seqinfo_ptr(n);
	seqs_processed += n;
}

}}