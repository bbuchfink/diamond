/****
DIAMOND protein aligner
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

#include <inttypes.h>
#include <cstdio>
#include <memory>
#include "multinode.h"
#include "data/block/block.h"
#include "search/lin_index/lin_index.h"
#include "search/search.h"
#include "util/log_stream.h"
#include "util/parallel/multiprocessing.h"
#include "data/fasta/volumed_fasta_file.h"

using std::runtime_error;
using std::string;
using std::unique_ptr;
using std::vector;

bool use_lin_index(const Job& job) {
	return job.is_linear_round() && !config.mutual_cover.present();
}

string lin_index_file(const string& volume_path) {
	return volume_path + ".seedidx";
}

void build_lin_indices(Job& job, const VolumedFile& volumes) {
	log_rss();
	const string base_dir = job.base_dir() + PATH_SEPARATOR + "seed_index" + PATH_SEPARATOR;
	if (config.recover) {
		recover_queue(job, "Seed index building", base_dir + "queue", base_dir + "finished", (int64_t)volumes.size());
		return;
	}
	job.make_temp_dir(base_dir);
	AtomicList finished(base_dir + "finished", job);
	Atomic q(base_dir + "queue", job, finished);
	const int64_t n = (int64_t)volumes.size();
	int64_t v, volumed_processed = 0;

	unique_ptr<vector<BitVector>> no_seed_hits;
	Search::Config cfg(no_seed_hits);
	Search::setup_search(config.sensitivity, cfg);

	/* Estimated number of pivots per block, written by the superblock partitioning. The
	   hash tables of the index are sized from it, which is what saves the pass over the
	   seeds of the block that would otherwise be needed to size them. */
	const auto seed_counts = read_seed_counts(seed_count_file(volumes.list_file()));

	while (job.goon() && (v = q.fetch_add(), v < n)) {
		const string index_file = lin_index_file(volumes[v].path);
		if (file_exists(index_file)) {
			finished.push(v);
			continue;
		}
		job.log("Building seed index. Block=%lli/%lli", v + 1, n);
		// The positions stored in the index refer to the order of the block on disk,
		// which is left as it is: the index keeps the length of the sequence a seed
		// comes from, so the pivot of a seed can be selected without sorting the block
		// by decreasing length first.
		TaskTimer timer("Loading block");
		unique_ptr<SequenceFile> file(new VolumedFastaFile({ volumes[v].path }, SequenceFile::Flags::SEQS, amino_acid_traits));
		unique_ptr<Block> block(file->load_seqs(INT64_MAX));
		file->close();
		timer.finish();
		// The index is written under a temporary name and renamed, so that a
		// worker crashing mid-write does not leave a truncated index behind.
		const string tmp_file = index_file + "." + std::to_string(job.worker_id()) + ".tmp";
		const auto counts = seed_counts.find(volumes[v].path);
		if (counts == seed_counts.end())
			throw runtime_error("No seed counter found for block " + volumes[v].path);
		Search::build_lin_index(*block, tmp_file, cfg, config.threads_, counts->second);
		block.reset();
		log_rss();
		if (!file_exists(index_file) && std::rename(tmp_file.c_str(), index_file.c_str()) != 0)
			throw runtime_error("Error renaming seed index file " + tmp_file);
		std::remove(tmp_file.c_str());
		finished.push(v);
		++volumed_processed;
	}
	if (volumed_processed > 0)
		job.finish_step();
	if (!job.goon())
		return;
	if (!finished.await(n, job))
		return;
	rmdir(base_dir.c_str());
}

void remove_lin_indices(const VolumedFile& volumes) {
	for (const Volume& v : volumes)
		remove_tmp_file(lin_index_file(v.path));
	remove_tmp_file(volumes.dir() + PATH_SEPARATOR + "seed_counts.tsv");
}