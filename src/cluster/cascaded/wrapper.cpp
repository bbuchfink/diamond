/****
DIAMOND protein aligner
Copyright (C) 2019-2023 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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

#include <numeric>
#include "cascaded.h"
#include "../data/fasta/fasta_file.h"
#include "../run/workflow.h"
#include "../util/system/system.h"
#include "../../search/search.h"

using std::shared_ptr;
using std::endl;
using std::tie;
using std::pair;
using std::iota;
using std::string;
using std::vector;
using std::ostream;
using std::ofstream;
using std::numeric_limits;
using std::tuple;
using std::get;
using std::ignore;
using std::back_inserter;
using std::unique_ptr;
using std::bind;
using std::function;
using namespace Util::Tsv;

namespace Cluster {
		
struct Config {
	Config(shared_ptr<SequenceFile>& db) :
		linclust(config.command == ::Config::LINCLUST),
		message_stream(true),
		verbosity(1),
		sens(from_string<Sensitivity>(rstrip(cluster_steps(config.approx_min_id, linclust).back(), "_lin"))),
		output_format(init_output(-1)),
		centroids(new FastaFile("", true, FastaFile::WriteAccess())),
		seqs_processed(0),
		letters_processed(0),
		oid_to_centroid_oid(new File(Schema{ Type::INT64, Type::INT64 }, "", Flags::TEMP))
	{
	}
	bool                                linclust;
	MessageStream                       message_stream;
	int                                 verbosity;
	Sensitivity                         sens;
	double                              block_size;
	std::unique_ptr<OutputFormat>       output_format;
	std::shared_ptr<FastaFile>          centroids;
	TaskTimer                          total_time;
	int64_t                             seqs_processed;
	int64_t                             letters_processed;
	std::vector<OId>                    centroid2oid;
	std::unique_ptr<File>               oid_to_centroid_oid;
};

int64_t seq_mem_use(Loc len, Loc id_len, int c, int min) {
	assert(min > 1);
	int extend_stage = (len / (min / 2)) * (15  // trace point buffer
		+ 16)        // SeedHitList
		+ 12         // SeedHitList
		+ 2 * len;   // SequenceSet
	extend_stage /= 2;
	assert(min > 1);
	return std::max(
		len + 8  // SequenceSet
		+ id_len + 8 // Seqtitle
		+ len * 9 / c / (min / 2) // SeedArray
		+ 8 // super_block_id_to_oid
		+ 8 // BestCentroid / clustering
		+ 4 // unaligned 
		, extend_stage);
}

struct BestCentroid : public Consumer, public vector<OId> {
	BestCentroid(OId block_size) :
		vector<OId>(block_size, -1)
	{}
	virtual void consume(const char *ptr, size_t n) {
		const char* const end = ptr + n;
		while (ptr < end) {
			const auto edge = *(Output::Format::Edge::Data*)ptr;
			ptr += sizeof(Output::Format::Edge::Data);
			this->operator[](edge.query) = edge.target;
		}
	}
	virtual void finalize() {}
	virtual ~BestCentroid() = default;
};

static vector<SuperBlockId> search_vs_centroids(shared_ptr<FastaFile>& super_block, const OId* super_block_id_to_oid, Config& cfg) {
	//if (cfg.verbosity >= 2)
		message_stream << "Searching vs. centroids #sequences = " << super_block->sequence_count() << " , #centroids = " << cfg.centroids->sequence_count() << endl;
	
	config.output_format = { "edge" };
	config.self = false;
	config.max_target_seqs_ = 1;
	config.toppercent = 100;
	config.sensitivity = cfg.sens;
	config.query_cover = config.member_cover;
	config.subject_cover = 0;
	config.query_or_target_cover = 0;
	if (cfg.linclust) {
		config.iterate.unset();
		config.linsearch = true;
	}
	else {
		config.iterate = vector<string>();
		config.linsearch = false;
	}
	config.lin_stage1 = false;
	tie(config.chunk_size, config.lowmem_) = block_size(Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT)), cfg.sens, cfg.linclust);
	cfg.centroids->set_seqinfo_ptr(0);
	shared_ptr<BestCentroid> best_centroid(new BestCentroid(super_block->sequence_count()));
	log_rss();
	Search::run(cfg.centroids, super_block, best_centroid);

	int64_t clustered = 0;
	vector<SuperBlockId> unaligned;
	for (SuperBlockId i = 0; i < super_block->sequence_count(); ++i) {
		const OId oid = super_block_id_to_oid[i];
		if (best_centroid->operator[](i) == -1) {
			unaligned.push_back(i);
		}
		else {
			const OId centroid_oid = cfg.centroid2oid[best_centroid->operator[](i)];
			cfg.oid_to_centroid_oid->write_record(centroid_oid, oid);
			++clustered;
		}
	}

	//if (cfg.verbosity >= 2)
		cfg.message_stream << clustered << " sequences assigned to clusters, " << unaligned.size() << " unaligned." << endl;
	return unaligned;
}

void Cascaded::run() {
	config.database.require();
	init_thresholds();
	config.hamming_ext = config.approx_min_id >= 50.0;
	TaskTimer total_time;
	TaskTimer timer("Opening the input file");
	shared_ptr<SequenceFile> db(SequenceFile::auto_create({ config.database }, SequenceFile::Flags::NEED_LETTER_COUNT | SequenceFile::Flags::OID_TO_ACC_MAPPING));
	if (db->type() == SequenceFile::Type::BLAST)
		throw std::runtime_error("Clustering is not supported for BLAST databases.");
	timer.finish();
	message_stream << "Input database: " << db->file_name() << " (" << db->sequence_count() << " sequences, " << db->letters() << " letters)" << endl;
	const int64_t mem_limit = Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT));
	const int64_t block_size = (int64_t)(::block_size(mem_limit, Sensitivity::FASTER, true).first * 1e9);
	unique_ptr<Util::Tsv::File> out(open_out_tsv());

	if (block_size >= (double)db->letters() && db->sequence_count() < numeric_limits<SuperBlockId>::max()) {
		const auto centroids = cascaded(db, config.command == ::Config::LINCLUST);
		timer.go("Generating output");
		output_mem<SuperBlockId>(*out, *db, centroids);
	}
	else {
		timer.go("Length sorting the input file");
		Config cfg(db);
		config.db_size = db->letters();
		const int minimizer_window = std::max(Search::sensitivity_traits[(int)align_mode.sequence_type].at(Sensitivity::FASTER).minimizer_window, 1);
		auto seq_size = function<int64_t(Loc)>(bind(seq_mem_use, std::placeholders::_1, 0, 1, minimizer_window));
		vector<tuple<FastaFile*, vector<OId>, Util::Tsv::File*>> super_blocks = db->length_sort(mem_limit / 2, seq_size);
		timer.finish();
		config.freq_masking = true;
		int i = 0;
		for (auto& b : super_blocks) {
			message_stream << "Processing super block " << (i++) + 1 << '/' << super_blocks.size() << endl;
			log_rss();
			log_stream << "Mem_sizes " << db->mem_size() << ' ' << cfg.centroids->mem_size() << endl;
			vector<OId> super_block_id_to_oid;
			Util::Tsv::File* oid_mapping;
			tie(ignore, super_block_id_to_oid, oid_mapping) = b;
			shared_ptr<FastaFile> seqs(get<0>(b));
			shared_ptr<SequenceFile> unaligned_db;
			vector<SuperBlockId> unaligned;
			timer.go("Reading super block mapping file");
			oid_mapping->template read<int64_t>(back_inserter(super_block_id_to_oid));
			timer.finish();
			log_rss();
			if (i == 1) {
				unaligned_db = seqs;
				unaligned.resize(seqs->sequence_count());
				iota(unaligned.begin(), unaligned.end(), 0);
			}
			else {
				unaligned = search_vs_centroids(seqs, super_block_id_to_oid.data(), cfg);
				timer.go("Creating subdatabase");
				unaligned_db.reset(seqs->sub_db(unaligned.cbegin(), unaligned.cend()));
				timer.go("Freeing memory");
				seqs->close();
				seqs.reset();
				timer.finish();
			}
			const vector<SuperBlockId> clustering = cascaded(unaligned_db, cfg.linclust);
			timer.go("Updating clustering");
			vector<SuperBlockId> centroids;
			for (SuperBlockId i = 0; i < (SuperBlockId)unaligned.size(); ++i) {
				const OId member_oid = super_block_id_to_oid[unaligned[i]], centroid_oid = super_block_id_to_oid[unaligned[clustering[i]]];
				cfg.oid_to_centroid_oid->write_record(centroid_oid, member_oid);
				if (member_oid == centroid_oid) {
					cfg.centroid2oid.push_back(centroid_oid);
					centroids.push_back(i);
				}
			}
			unaligned_db->sub_db(centroids.cbegin(), centroids.cend(), cfg.centroids.get());
			timer.go("Freeing memory");
			unaligned_db->close();
			unaligned_db.reset();
			delete oid_mapping;
			timer.finish();
		}
		message_stream << "Total clusters: " << cfg.centroid2oid.size() << endl;
		message_stream << "Total time: " << total_time.seconds() << 's' << endl;
		timer.go("Generating output");
		output_mem(*out, *db, *cfg.oid_to_centroid_oid);
	}
	db->close();
}

}