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

#include "util/memory/memory_resource.h"
#include "search/search.h"
#include "external.h"
#include "file_array.h"
#include "data/sequence_file.h"
#include "basic/shape.h"
#include "basic/seed_iterator.h"
#include "basic/shape_config.h"
#include "build_pair_table.h"
#include "util/string/string.h"
#include "masking/masking.h"
#include "search/seed_complexity.h"

using std::pmr::monotonic_buffer_resource;
using std::atomic;
using std::unique_ptr;
using std::vector;
using std::string;
using std::ofstream;
using std::runtime_error;

namespace External {

RadixedTable build_seed_table(Job& job, const VolumedFile& volumes, int shape) {
	const int64_t BUF_SIZE = 4096;
	const double SKETCH_SIZE_RATIO = 0.1;
	MaskingAlgo masking = MaskingAlgo::NONE;
	if (config.soft_masking.empty() || config.soft_masking == "tantan")
		masking |= MaskingAlgo::TANTAN;
	if(config.motif_masking.empty() || config.motif_masking == "1")
		masking |= MaskingAlgo::MOTIF;
	//Reduction::set_reduction(Search::no_reduction);
	const auto traits = Search::sensitivity_traits.at(config.sensitivity);
	Loc sketch_size = config.sketch_size == 0 ? traits.sketch_size : config.sketch_size;
	const double seed_cut = config.seed_cut_ == 0.0 ? traits.seed_cut : config.seed_cut_,
		seed_complexity_cut = seed_cut * std::log(2.0) * ::shapes[0].weight_;
	if (sketch_size == 0)
		sketch_size = std::numeric_limits<Loc>::max();
	//const uint64_t shift = shapes[shape].bit_length() - RADIX_BITS;

	const std::string base_dir = job.base_dir() + PATH_SEPARATOR + "seed_table_" + std::to_string(shape) + PATH_SEPARATOR, qpath = base_dir + "queue";
	if (job.round() == 0 && shape == 0)
		mkdir(job.root_dir() + "accessions");
	mkdir(base_dir);
	unique_ptr<FileArray> output_files(new FileArray(base_dir, RADIX_COUNT, job.worker_id(), false));

	Atomic q(qpath, job);
	atomic<int> volumes_processed(0);
	SimpleThreadPool pool;
	ClusterStats stats_all;
	auto worker = [&](const atomic<bool>& stop, int thread_id) {
		const bool accs = job.round() == 0 && shape == 0;
		BufferArray buffers(*output_files, RADIX_COUNT);
		monotonic_buffer_resource mem_pool;
		int64_t v = 0;
		ClusterStats stats;
		while (!stop.load(std::memory_order_relaxed) && (v = q.fetch_add(), v < (int64_t)volumes.size())) {
			job.log("Building seed table. Shape=%i/%i Volume=%lli/%lli Records=%s", shape + 1, ::shapes.count(), v + 1, volumes.size(), Util::String::format(volumes[v].record_count).c_str());
			unique_ptr<SequenceFile> in(SequenceFile::auto_create({ volumes[v].path }));
			unique_ptr<ofstream> acc_out;
			if (accs) {
				const string name = job.root_dir() + "accessions" + PATH_SEPARATOR + std::to_string(v) + ".txt";
				acc_out.reset(new ofstream(name));
				if (!acc_out->good())
					throw runtime_error("Error opening file " + name);
			}
			string id;
			vector<Letter> seq;
			int64_t oid = volumes[v].oid_begin;
			while (!stop.load(std::memory_order_relaxed) && in->read_seq(seq, id, nullptr)) {
				if (job.round() > 0)
					oid = atoll(id.c_str());
				if (accs)
					*acc_out << Util::Seq::seqid(id.c_str()) << std::endl;
				if (masking != MaskingAlgo::NONE) {
					MaskingTable table;
					stats.masking_stat += Masking::get()(seq.data(), seq.size(), masking, 0, &table);
				}
				std::pmr::vector<Letter> buf = Reduction::reduce_seq(Sequence(seq), mem_pool);
				const Shape& sh = shapes[shape];
				if (seq.size() < (size_t)sh.length_) {
					++oid;
					continue;
				}
				//SketchIterator it(buf, sh, std::min(sketch_size, std::max((int)std::round(seq.size() * SKETCH_SIZE_RATIO), 1)));
				SketchIterator it(buf.cbegin(), buf.cend(), sh, sketch_size);
				while (it.good()) {
					++stats.seeds_considered;
					if(Search::seed_is_complex(&seq[it.pos()], sh, seed_complexity_cut)) {
						const uint64_t key = *it;
						buffers.write_msb(SeedEntry(key, oid, (int32_t)seq.size()));
						++stats.seeds_indexed;
					}
					++it;
				}
				++oid;
			}
			in->close();
			volumes_processed.fetch_add(1, std::memory_order_relaxed);
		}
		stats_all.add(stats);
		};
	for (int i = 0; i < config.threads_; ++i)
		pool.spawn(worker, i);
	pool.join_all();
	const RadixedTable buckets = output_files->buckets();
	TaskTimer timer("Closing the output files");
	output_files.reset();
	Atomic finished(base_dir + "finished", job);
	finished.fetch_add(volumes_processed);
	finished.await(volumes.size());
	job.log(stats_all);
	return buckets;
}

}