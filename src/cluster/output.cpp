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

#include <mutex>
#include <memory>
#include "cluster.h"
#include "util/parallel/thread_pool.h"
#include "util/algo/merge_files.h"
#include "dp/dp.h"
#include "stats/hauser_correction.h"
#include "output/output.h"
#include "align/culling.h"

using std::vector;
using std::unique_ptr;
using std::endl;
using std::to_string;
using std::mutex;
using std::lock_guard;
using std::thread;
using std::make_pair;
using std::list;
using std::lower_bound;
using std::atomic;
using std::shared_ptr;
using std::function;
using std::string;

namespace Cluster {
	
struct OutputWriter {
	OutputWriter(File* file_, char sep = char(0), bool first = true) :
		file_(file_),
		first(first),
		sep(sep)
	{
	};
	void operator()(TextBuffer* buf) {
		if (!first && sep != '\0')
		{
			file_->write(sep);
		}
		file_->write(buf->data(), buf->size());
		first = false;
	}
	File* file_;
	bool first;
	char sep;
};

static void align_centroid(CentroidId centroid, ReorderQueue<TextBuffer*, Cluster::OutputWriter>& out, Statistics& stats, ThreadPool& tp, Cfg& cfg) {
	DP::Targets dp_targets;
	const OId centroid_oid = cfg.centroids[centroid];
	const BlockId centroid_id = cfg.centroid_block->oid2block_id(centroid_oid);
	const Sequence centroid_seq = cfg.centroid_block->seqs()[centroid_id];
	auto begin = lower_bound(cfg.clusters.cbegin(centroid), cfg.clusters.cend(centroid), cfg.member_block->oid_begin());
	auto end = lower_bound(cfg.clusters.cbegin(centroid), cfg.clusters.cend(centroid), cfg.member_block->oid_end());

	for (auto it = begin; it != end; ++it) {
		const BlockId block_id = cfg.member_block->oid2block_id(*it);
		const Sequence seq(cfg.member_block->seqs()[block_id]);
		const int bin = DP::BandedSwipe::bin(cfg.hsp_values, centroid_seq.length(), 0, 0, (int64_t)seq.length() * (int64_t)centroid_seq.length(), 0, 0);
		dp_targets[bin].emplace_back(seq, seq.length(), block_id);
	}

	const HauserCorrection cbs(centroid_seq);
	const string centroid_seqid = cfg.lazy_titles ? cfg.db->seqid(centroid_oid, false, false) : cfg.centroid_block->ids()[centroid_id];
	DP::Params p{ centroid_seq, centroid_seqid.c_str(), Frame(0), centroid_seq.length(), config.comp_based_stats == 1 ? cbs.int8.data() : nullptr, DP::Flags::FULL_MATRIX, false, 0, 0, cfg.hsp_values, stats, &tp };
	list<Hsp> hsps = DP::BandedSwipe::swipe(dp_targets, p);

	TextBuffer* buf = new TextBuffer;
	const int replen = centroid_seq.length();
	hsps.sort(Hsp::OrderBySwipeTarget());
	for (Hsp& hsp : hsps) {
		const Loc memberlen = cfg.member_block->seqs().length(hsp.swipe_target);
		if (Extension::filter_hsp(hsp, replen, nullptr, memberlen, nullptr, centroid_seq, Sequence(), 0, 0, nullptr))
			continue;
		serialize(HspContext(hsp,
			centroid_id,
			centroid_oid,
			TranslatedSequence(centroid_seq),
			centroid_seqid.c_str(),
			cfg.member_block->block_id2oid(hsp.swipe_target),
			memberlen,
			cfg.lazy_titles ? cfg.db->seqid(cfg.member_block->block_id2oid(hsp.swipe_target), false, false).c_str() : cfg.member_block->ids()[hsp.swipe_target],
			0,
			0,
			Sequence()), *buf);
	}
	out.push(centroid, buf);
}

File* realign_block_pair(CentroidId begin, CentroidId end, Cfg& cfg) {
	File* out = new File { Temporary() };
	Cluster::OutputWriter writer(out);
	ReorderQueue<TextBuffer*, Cluster::OutputWriter> output_sink(begin, writer);
	auto worker = [&](ThreadPool& tp, int64_t i) {
		Statistics stats;
		align_centroid(i, output_sink, stats, tp, cfg);
		statistics += stats;
	};
	ThreadPool tp(worker, begin, end);
	tp.run(config.threads_, true);
	tp.join();
	out->rewind();
	return out;
}

void realign(const FlatArray<OId>& clusters, const vector<OId>& centroids, SequenceFile& db, function<void(const HspContext&)>& callback, HspValues hsp_values) {
	const int64_t block_size = Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT)) / 2;
	*message_stream << "Block size: " << block_size << " byte." << endl;
	db.set_seqinfo_ptr(0);
	score_matrix.set_db_letters(config.db_size ? config.db_size : db.letters().value());
	OId centroid_offset = 0;
	int n1 = 0;
	TaskTimer timer;
	Cfg cfg{ hsp_values, flag_any(db.format_flags(), SequenceFile::FormatFlags::TITLES_LAZY), clusters, centroids, &db };

	db.flags() |= SequenceFile::Flags::SEQS;
	if (!cfg.lazy_titles)
		db.flags() |= SequenceFile::Flags::TITLES;

	while (centroid_offset < db.sequence_count().value()) {
		timer.go("Loading centroid block");
		db.set_seqinfo_ptr(centroid_offset);
		cfg.centroid_block.reset(db.load_seqs(block_size, 0, nullptr));
		centroid_offset = db.tell_seq();

		db.set_seqinfo_ptr(0);
		int n2 = 0;
		const CentroidId begin = lower_bound(centroids.cbegin(), centroids.cend(), cfg.centroid_block->oid_begin()) - centroids.cbegin(),
			end = lower_bound(centroids.cbegin(), centroids.cend(), cfg.centroid_block->oid_end()) - centroids.cbegin();
		timer.finish();
		*log_stream << "Total centroids = " << end - begin << endl;
		vector<File*> tmp;
		for (;;) {
			if (cfg.centroid_block->seqs().size() == db.sequence_count().value())
				cfg.member_block = cfg.centroid_block;
			else {
				timer.go("Loading member block");
				cfg.member_block.reset(db.load_seqs(block_size, 0, nullptr));
			}
			if (cfg.member_block->empty())
				break;
			timer.go("Processing centroid block " + to_string(n1 + 1) + ", member block " + to_string(n2 + 1));
			tmp.push_back(realign_block_pair(begin, end, cfg));
			++n2;
			if (cfg.centroid_block->seqs().size() == db.sequence_count().value())
				break;
		}
		timer.go("Joining centroid block " + to_string(n1 + 1));
	
		merge_sorted_files<HspContext, vector<File*>::iterator, decltype(callback)>(tmp.begin(), tmp.end(), callback);
		for (File* f : tmp) {
			f->close();
			delete f;
		}
		++n1;
	}
	timer.finish();
	statistics.print();
}

void realign(const vector<OId>& clustering, SequenceFile& db, function<void(const HspContext&)>& callback, HspValues hsp_values) {
	TaskTimer timer("Finding clusters");
	FlatArray<OId> clusters;
	vector<OId> centroids;
	tie(clusters, centroids) = cluster_sorted(clustering);
	timer.finish();
	realign(clusters, centroids, db, callback, hsp_values);
}

}