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

#include <inttypes.h>
#include <fstream>
#include <unordered_map>
#include "multinode.h"
#include "volume.h"
#include "basic/config.h"
#include "util/memory/memory_resource.h"
#include "util/io/compressed_buffer.h"
#include "file_array.h"
#include "input_buffer.h"

using std::endl;
using std::unordered_map;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::runtime_error;
using std::unique_ptr;

static OId output_oids(Job& job, const vector<OId>& merged) {
	FILE* out = fopen(config.output_file.c_str(), "wt");
	if (!out)
		throw runtime_error("Error opening output file: " + config.output_file);
	OId n = 0;
	for (OId i = 0; i <= job.max_oid(); ++i) {
		if (merged[i] == i) ++n;
		fprintf(out, "%" PRId64 "\t%" PRId64 "\n", (int64_t)merged[i], (int64_t)i);
	}
	fclose(out);
	return n;
}

struct AccMapping {
	static constexpr bool POD = false;
	static constexpr OId NIL = std::numeric_limits<OId>::max();
	OId rep = NIL, member = NIL;
	std::pmr::string rep_acc, member_acc;
	AccMapping(OId rep, OId member, std::pmr::string&& member_acc, std::pmr::memory_resource& pool) :
		rep(rep),
		member(member),
		rep_acc(&pool),
		member_acc(std::move(member_acc))
	{
	}
	AccMapping(OId rep, OId member, const std::pmr::string& member_acc, std::pmr::memory_resource& pool) :
		rep(rep),
		member(member),
		rep_acc(&pool),
		member_acc(member_acc, &pool)
	{
	}
	AccMapping(std::pmr::memory_resource& pool) :
		rep(NIL),
		member(NIL),
		rep_acc(&pool),
		member_acc(&pool)
	{
	}
	OId key() const {
		return rep;
	}
	bool operator<(const AccMapping& m) const {
		return rep < m.rep || (rep == m.rep && member < m.member);
	}
	friend void serialize(const AccMapping& m, CompressedBuffer& buf) {
		buf.write(m.rep);
		buf.write(m.member);
		buf.write(m.rep_acc.c_str(), m.rep_acc.length() + 1);
		buf.write(m.member_acc.c_str(), m.member_acc.length() + 1);
	}
	friend void deserialize(InputFile& f, AccMapping& m) {
		f.read(&m.rep);
		f.read(&m.member);
		f >> m.rep_acc;
		f >> m.member_acc;
	}
};

static std::pmr::unordered_map<OId, std::pmr::string> read_mapping_table(Job& job, const Volume& vol, size_t v, std::pmr::memory_resource& pool, bool remove) {
	std::pmr::unordered_map<OId, std::pmr::string> oid2acc(&pool);
	oid2acc.reserve(vol.record_count);
	const string path = job.root_dir() + "input" + std::to_string(v) + ".tsv";
	ifstream in(path);
	if (!in.good())
		throw runtime_error("Error opening accessions file: " + path);
	OId oid;
	std::pmr::string acc(&pool);
	while (in >> oid) {
		in >> acc;
		if (!in)
			throw runtime_error("Format error in accessions file: " + path);
		if (oid2acc.emplace(oid, acc).second == false)
			throw runtime_error("Duplicate OID in accessions file: " + path);
	}
	in.close();
	if (oid2acc.size() != vol.record_count)
		throw runtime_error("Accessions file does not contain all OIDs");
	if(remove)
		remove_tmp_file(path);
	return oid2acc;
}

static RadixedTable output_round1(Job& job, const vector<OId>& merged, const VolumedFile& volumes) {
	std::pmr::unsynchronized_pool_resource pool;
	const string base_dir = job.root_dir() + "output" + PATH_SEPARATOR;
	job.make_temp_dir(base_dir);
	unique_ptr<FileArray> output_files(new FileArray(base_dir, RADIX_COUNT, job.worker_id(), false));
	const int shift = std::max(bit_length(job.max_oid()) - RADIX_BITS, 0);
	BufferArray buffers(*output_files, RADIX_COUNT);
	for (size_t v = 0; v < volumes.size(); ++v) {
		const Volume& vol = volumes.at(v);
		job.log("Building output table (round 1) volume %zu/%zu", v + 1, volumes.size());
		const std::pmr::unordered_map<OId, std::pmr::string> oid2acc = read_mapping_table(job, vol, v, pool, false);
		for (auto it = oid2acc.cbegin(); it != oid2acc.cend(); ++it) {
			const OId oid = it->first;
			const std::pmr::string& acc = it->second;
			AccMapping m(merged[oid], oid, acc, pool);
			buffers.write(m.rep >> shift, m);
		}
		log_rss();
	}
	return output_files->buckets(shift);
}

static OId output_round2(Job& job, const vector<OId>& merged, const VolumedFile& volumes, const RadixedTable& round1) {
	std::pmr::unsynchronized_pool_resource pool;
	ofstream out(config.output_file);
	if (!out.good())
		throw runtime_error("Error opening output file: " + config.output_file);
	OId cluster_count = 0;
	for (size_t v = 0; v < volumes.size(); ++v) {
		const Volume& vol = volumes.at(v);
		job.log("Building output table (round 2) volume %zu/%zu oid %" PRId64 " - %" PRId64, v + 1, volumes.size(), vol.oid_begin, vol.oid_end);
		const std::pmr::unordered_map<OId, std::pmr::string> oid2acc = read_mapping_table(job, vol, v, pool, true);
		for (const Bucket& b : round1) {
			log_rss();
			if (b.key_end() <= vol.oid_begin || b.key_begin() >= vol.oid_end)
				continue;
			VolumedFile f(b);
			InputBuffer<AccMapping> data(f, pool);
			data.sort();
			OId last_rep = AccMapping::NIL;
			for (auto it = data.cbegin(); it != data.cend(); ++it) {
				if (it->rep < vol.oid_begin || it->rep >= vol.oid_end)
					continue;
				out << oid2acc.at(it->rep) << '\t' << it->member_acc << endl;
				if(!out.good())
					throw runtime_error("Error writing output file: " + config.output_file);
				if (it->rep != last_rep) {
					++cluster_count;
					last_rep = it->rep;
				}
			}
		}
	}
	for (const Bucket& b : round1)
		VolumedFile(b).remove();
	rmdir(job.root_dir() + "output");
	return cluster_count;
}

void merge(Job& job, const VolumedFile& volumes) {
	job.log("Merging clusterings");
	const vector<OId> merged = build_merged(job);
	OId n;
	if (config.oid_output)
		n = output_oids(job, merged);
	else {
		const RadixedTable round1 = output_round1(job, merged, volumes);
		n = output_round2(job, merged, volumes, round1);
	}
	job.log("Total clusters: %" PRId64, (int64_t)n);
}