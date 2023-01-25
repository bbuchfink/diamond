/****
DIAMOND protein aligner
Copyright (C) 2021-2022 Max Planck Society for the Advancement of Science e.V.

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

#include <mutex>
#include <memory>
#include "cluster.h"
#include "../util/algo/algo.h"
#include "../util/util.h"
#include "../util/algo/sort_helper.h"
#include "../util/parallel/thread_pool.h"
#include "../dp/dp.h"
#include "../output/output_format.h"
#include "../stats/hauser_correction.h"
#include "../output/output.h"
#include "../util/algo/merge_files.h"

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

struct Cfg {
	Cfg(HspValues hsp_values, bool lazy_titles, const FlatArray<OId>& clusters, const vector<OId>& centroids, SequenceFile& db):
		hsp_values(hsp_values),
		lazy_titles(lazy_titles),
		clusters(clusters),
		centroids(centroids),		
		db(db)
	{}
	const HspValues hsp_values;
	const bool lazy_titles;
	const FlatArray<OId>& clusters;
	const vector<OId>& centroids;
	SequenceFile& db;
	shared_ptr<Block> centroid_block, member_block;
};

static void align_centroid(CentroidId centroid, ReorderQueue<TextBuffer*, OutputWriter>& out, Statistics& stats, ThreadPool& tp, Cfg& cfg) {
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

	const Bias_correction cbs(centroid_seq);
	const string centroid_seqid = cfg.lazy_titles ? cfg.db.seqid(centroid_oid) : cfg.centroid_block->ids()[centroid_id];
	DP::Params p{ centroid_seq, centroid_seqid.c_str(), Frame(0), centroid_seq.length(), config.comp_based_stats == 1 ? cbs.int8.data() : nullptr, DP::Flags::FULL_MATRIX, cfg.hsp_values, stats, &tp };
	list<Hsp> hsps = DP::BandedSwipe::swipe(dp_targets, p);

	TextBuffer* buf = new TextBuffer;
	TypeSerializer<HspContext> s(*buf);
	for (Hsp& hsp : hsps) {
		s << HspContext(hsp,
			centroid_id,
			centroid_oid,
			TranslatedSequence(centroid_seq),
			centroid_seqid.c_str(),
			cfg.member_block->block_id2oid(hsp.swipe_target),
			cfg.member_block->seqs().length(hsp.swipe_target),
			cfg.lazy_titles ? cfg.db.seqid(cfg.member_block->block_id2oid(hsp.swipe_target)).c_str() : cfg.member_block->ids()[hsp.swipe_target],
			0,
			0,
			Sequence());
	}
	out.push(centroid, buf);
}

static InputFile* run_block_pair(CentroidId begin, Cfg& cfg) {
	atomic<CentroidId> next(begin);
	TempFile out;
	OutputWriter writer{ &out };
	output_sink.reset(new ReorderQueue<TextBuffer*, OutputWriter>(begin, writer));
	auto worker = [&](ThreadPool& tp) {
		Statistics stats;
		CentroidId i = next++;
		if (i >= (CentroidId)cfg.centroids.size() || cfg.centroids[i] >= cfg.centroid_block->oid_end())
			return false;
		align_centroid(i, *output_sink, stats, tp, cfg);
		statistics += stats;
		return true;
	};
	ThreadPool tp(worker);
	tp.run(config.threads_, true);
	tp.join();
	InputFile* f = new InputFile(out);
	return f;
}

void realign(const FlatArray<OId>& clusters, const vector<OId>& centroids, SequenceFile& db, function<void(const HspContext&)>& callback, HspValues hsp_values) {
	const int64_t block_size = Util::String::interpret_number(config.memory_limit.get(DEFAULT_MEMORY_LIMIT)) / 2;
	message_stream << "Block size: " << block_size << " byte." << endl;
	db.set_seqinfo_ptr(0);
	score_matrix.set_db_letters(config.db_size ? config.db_size : db.letters());
	OId centroid_offset = 0;
	int n1 = 0;
	task_timer timer;
	Cfg cfg{ hsp_values, flag_any(db.format_flags(), SequenceFile::FormatFlags::TITLES_LAZY), clusters, centroids, db };

	SequenceFile::LoadFlags flags = SequenceFile::LoadFlags::SEQS | SequenceFile::LoadFlags::CONVERT_ALPHABET | SequenceFile::LoadFlags::NO_CLOSE_WEAKLY;
	if (!cfg.lazy_titles)
		flags |= SequenceFile::LoadFlags::TITLES;

	while (centroid_offset < db.sequence_count()) {
		timer.go("Loading centroid block");
		db.set_seqinfo_ptr(centroid_offset);
		cfg.centroid_block.reset(db.load_seqs(block_size, nullptr, flags));
		centroid_offset = db.tell_seq();

		db.set_seqinfo_ptr(0);
		int n2 = 0;
		const CentroidId begin = lower_bound(centroids.cbegin(), centroids.cend(), cfg.centroid_block->oid_begin()) - centroids.cbegin(),
			end = lower_bound(centroids.cbegin(), centroids.cend(), cfg.centroid_block->oid_end()) - centroids.cbegin();
		timer.finish();
		log_stream << "Total centroids = " << end - begin << endl;
		vector<InputFile*> tmp;
		for (;;) {
			if (cfg.centroid_block->seqs().size() == db.sequence_count())
				cfg.member_block = cfg.centroid_block;
			else {
				timer.go("Loading member block");
				cfg.member_block.reset(db.load_seqs(block_size, nullptr, flags));
			}
			if (cfg.member_block->empty())
				break;
			timer.go("Processing centroid block " + to_string(n1 + 1) + ", member block " + to_string(n2 + 1));
			tmp.push_back(run_block_pair(begin, cfg));
			++n2;
			if (cfg.centroid_block->seqs().size() == db.sequence_count())
				break;
		}
		timer.go("Joining centroid block " + to_string(n1 + 1));
	
		merge_sorted_files<HspContext, vector<InputFile*>::iterator, decltype(callback)>(tmp.begin(), tmp.end(), callback);
		for (InputFile* f : tmp) {
			f->close_and_delete();
			delete f;
		}
		++n1;
	}
	timer.finish();
	statistics.print();
}

void realign(const vector<OId>& clustering, SequenceFile& db, function<void(const HspContext&)>& callback, HspValues hsp_values) {
	task_timer timer("Finding clusters");
	FlatArray<OId> clusters;
	vector<OId> centroids;
	tie(clusters, centroids) = cluster_sorted(clustering);
	timer.finish();
	realign(clusters, centroids, db, callback, hsp_values);
}

}