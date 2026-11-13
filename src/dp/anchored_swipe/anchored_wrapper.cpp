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

#include <memory>
#include <cfloat>
#include <mutex>
#include <iterator>
#include "../dp.h"
#include "config.h"
#include "anchored.h"
#include "../long_profile/score_profile.h"
#include "util/simd/dispatch.h"
#include "stats/score_matrix.h"
#include "align/def.h"
#include "stats/stats.h"
#include "util/memory/mem_profile.h"
#include "basic/config.h"

using std::vector;
using std::unique_ptr;
using std::pair;
using std::sort;

namespace DP { namespace AnchoredSwipe { namespace DISPATCH_ARCH {

// Returns the DP target of an extension, without its query letters or profile, which are set once it is loaded into a
// channel of the kernel.
static DP::AnchoredSwipe::Target<int16_t> align_right(Loc query_len, Sequence target_seq, bool reverse, Loc i, Loc j, Loc d_begin, Loc d_end,
	int64_t target_idx,
	const DP::AnchoredSwipe::Config& cfg)
{
	// The diagonals [d_begin, d_end) are computed.
	Loc qlen = query_len - i, tlen = target_seq.length();
	const Loc d0 = Geo::clip_diag(Geo::diag_sub_matrix(d_begin, i, j), qlen, tlen),
		d1 = Geo::clip_diag(Geo::diag_sub_matrix(d_end - 1, i, j), qlen, tlen);
	tlen = std::min(tlen, Geo::j(qlen - 1, d0) + 1);
	assert(tlen > 0);
	assert(d1 >= d0);
	assert(d1 >= 0 && d1 + 1 - d0 <= cfg.band_max);

	const Sequence clipped_target = reverse ? target_seq.subseq(target_seq.length() - tlen, target_seq.length()) : target_seq.subseq(0, tlen);
	DP::AnchoredSwipe::Target<int16_t> t(clipped_target, d0, d1 + 1, i, qlen, target_idx, reverse);
	t.profile = t.profile_rev = nullptr;
	t.query = t.query_rev = nullptr;
	return t;
}

static DP::AnchoredSwipe::Target<int16_t> align_left(Loc query_len, Sequence target_seq, Loc i, Loc j, Loc d_begin, Loc d_end,
	int64_t target_idx,
	const DP::AnchoredSwipe::Config& cfg)
{
	const Loc qlen = query_len, tlen = target_seq.length();
	const Loc ir = qlen - 1 - i, jr = tlen - 1 - j;
	return align_right(qlen, target_seq.subseq(0, j + 1), true, ir, jr, Geo::rev_diag(d_end, qlen, tlen), Geo::rev_diag(d_begin, qlen, tlen), target_idx, cfg);
}

// Returns the DP target of one extension. Each extension yields exactly one target.
static DP::AnchoredSwipe::Target<int16_t> make_target(const ExtensionPipeline::Extension& e, int64_t target_idx, const DP::AnchoredSwipe::Config& cfg)
{
	const Loc qlen = e.query.length();
	if (e.direction == ExtensionPipeline::RIGHT)
		return align_right(qlen, e.target.subseq(e.j), false, e.i, e.j, e.d_min, e.d_max, target_idx, cfg);
	else
		return align_left(qlen, e.target, e.i, e.j, e.d_min, e.d_max, target_idx, cfg);
}

// Writes seq framed by `padding` padding letters on each side to out, reversed if reverse is set, and returns the pointer to its first letter.
static const Letter* pad_query(Sequence seq, bool reverse, int64_t padding, Letter* out) {
	std::fill(out, out + padding, DP::AnchoredSwipe::PADDING_LETTER);
	Letter* begin = out + padding;
	if (reverse)
		std::reverse_copy(seq.data(), seq.end(), begin);
	else
		std::copy(seq.data(), seq.end(), begin);
	std::fill(begin + seq.length(), begin + seq.length() + padding, DP::AnchoredSwipe::PADDING_LETTER);
	return begin;
}

/* Target source of the kernel (see TargetArray in anchored.h) feeding it the extensions of the queue. It first serves
   the batch of extensions popped from the queue by the caller, ordered by band. Once they are used up, every channel that
   becomes idle pops the next extension from the queue, so an extension submitted by a target task while the kernel
   is running is loaded as soon as a channel is free for it.
   This goes on until the queue has no extension for an idle channel. From then on the source does not ask the queue
   again, so the kernel completes the channels still loaded and returns. Having the kernel wait for further extensions
   instead would keep it running with part of its channels idle at the full cost per column, and when the queue
   runs dry it does so for all the kernels popping from it, which would then share the few extensions coming in. Giving up
   leaves these to the kernels that still have all their channels loaded, and frees the thread for the tasks submitting
   extensions. A new kernel is started once enough extensions have been queued (see ThreadPool).
   The padded query (or the profile) of an extension is built when it is loaded into a channel and lives in the
   storage of that channel until the channel is loaded again. An extension whose pair carries a score table (a
   compositionally adjusted matrix) is computed with it, the others with the unadjusted matrix.
   Computed extensions are passed to the callback in batches of FLUSH_SIZE, and the remainder when the kernel is done.
   Every batch is ordered by the sequence in which the extensions were popped from the queue. */
class ExtensionSource {
public:

	enum { FLUSH_SIZE = 64 };

	ExtensionSource(ExtensionPipeline::ExtensionQueue& queue, std::vector<ExtensionPipeline::Extension>&& batch, const DP::AnchoredSwipe::Config& cfg,
		const ExtensionPipeline::ExtensionCallback& callback) :
		queue_(queue),
		batch_(std::move(batch)),
		batch_pos_(0),
		shard_(ExtensionPipeline::ExtensionQueue::thread_shard()),
		popped_(0),
		cfg_(cfg),
		callback_(callback),
		padding_(DP::AnchoredSwipe::MATRIX_ROW_SCORES
			? std::max(DP::AnchoredSwipe::DISPATCH_ARCH::profile_padding<::DISPATCH_ARCH::ScoreVector<int16_t, 0>>(cfg.band_max), (int64_t)LongScoreProfile<int16_t>::DEFAULT_PADDING)
			: DP::AnchoredSwipe::DISPATCH_ARCH::profile_padding<::DISPATCH_ARCH::ScoreVector<int16_t, 0>>(cfg.band_max)),
		profile_time_(0),
		output_time_(0),
		record_bands_(!config.band_hist_file.empty())
	{
		TaskTimer timer;
		batch_targets_.reserve(batch_.size());
		for (size_t k = 0; k < batch_.size(); ++k)
			batch_targets_.push_back(make_target(batch_[k], (int64_t)k, cfg));
		cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE_ADD, timer.microseconds());

		timer.go();
		sort(batch_targets_.begin(), batch_targets_.end());
		cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE_SORT, timer.microseconds());

		if (DP::AnchoredSwipe::MATRIX_ROW_SCORES)
			DP::AnchoredSwipe::make_score_table(score_matrix.matrix8(), 32, score_table_);
	}

	~ExtensionSource() {
		cfg_.stats.inc(Statistics::TIME_PROFILE, profile_time_ / 1000);
		cfg_.stats.inc(Statistics::TIME_ANCHORED_SWIPE_OUTPUT, output_time_ / 1000);
		cfg_.stats.inc(Statistics::EXTENSIONS_REFILLED, popped_);
		if (record_bands_)
			ExtensionPipeline::add_band_histogram(band_hist_);
	}

	bool next(int channel, DP::AnchoredSwipe::Target<int16_t>& t) {
		if ((size_t)channel >= channels_.size())
			channels_.resize(channel + 1);
		Channel& c = channels_[channel];
		if (batch_pos_ < batch_targets_.size()) {
			t = batch_targets_[batch_pos_++];
			c.extension = batch_[t.target_idx];
			c.seq = t.target_idx;
		}
		else {
			if (!pop(c.extension))
				return false;
			c.seq = int64_t(batch_.size()) + popped_++;
			t = make_target(c.extension, c.seq, cfg_);
		}
		if (record_bands_)
			record_band(c.extension.diag_spread());
		load(c, t);
		return true;
	}

	void finish(int channel, int16_t score, Loc query_end, Loc target_end) {
		Channel& c = channels_[channel];
		c.extension.score = score;
		c.extension.query_len = query_end;
		c.extension.target_len = target_end;
		done_.emplace_back(c.seq, c.extension);
		if (done_.size() >= FLUSH_SIZE)
			flush();
	}

	// Pops the next extension of the queue to e. Returns false if none is available, and keeps doing so from then on
	// without asking the queue again (see above).
	bool pop(ExtensionPipeline::Extension& e) {
		if (queue_drained_)
			return false;
		if (queue_.pop(shard_, e))
			return true;
		queue_drained_ = true;
		return false;
	}

	// Passes the computed extensions to the callback.
	void flush() {
		if (done_.empty())
			return;
		TaskTimer timer;
		sort(done_.begin(), done_.end(), [](const Done& a, const Done& b) { return a.first < b.first; });
		out_.clear();
		out_.reserve(done_.size());
		for (const Done& d : done_)
			out_.push_back(d.second);
		done_.clear();
		callback_(out_);
		output_time_ += timer.nanoseconds();
	}

	DP::AnchoredSwipe::Options options() const {
		return DP::AnchoredSwipe::Options{ nullptr, nullptr, score_table_ };
	}

private:

	struct Channel {
		ExtensionPipeline::Extension extension;
		int64_t seq = 0;                        // position in the sequence the extensions were popped in
		vector<Letter> query;                   // MATRIX_ROW_SCORES mode: padded query of the extension
		LongScoreProfile<int16_t> profile;      // otherwise: its profile, reversed for LEFT extensions
	};

	using Done = pair<int64_t, ExtensionPipeline::Extension>;

	void record_band(Loc band) {
		if ((size_t)band >= band_hist_.size())
			band_hist_.resize((size_t)band + 1, 0);
		++band_hist_[band];
	}

	// Builds the padded query or profile of the extension of the channel, for the direction it extends into.
	void load(Channel& c, DP::AnchoredSwipe::Target<int16_t>& t) {
		MEM_SCOPE("dp/query-profile");
		TaskTimer timer;
		const ExtensionPipeline::Extension& e = c.extension;
		const bool reverse = e.direction == ExtensionPipeline::LEFT;
		if (DP::AnchoredSwipe::MATRIX_ROW_SCORES) {
			c.query.resize((size_t)e.query.length() + 2 * padding_);
			(reverse ? t.query_rev : t.query) = pad_query(e.query, reverse, padding_, c.query.data());
			t.score_table = e.pair->score_table;
		}
		else {
			assert(e.pair->score_table == nullptr); // adjusted matrices are only supported in MATRIX_ROW_SCORES mode
			c.profile = make_profile16(e.query, nullptr, padding_, &score_matrix);
			if (reverse)
				for (int l = 0; l < AMINO_ACID_COUNT; ++l)
					std::reverse(c.profile.data[l].begin(), c.profile.data[l].end());
			(reverse ? t.profile_rev : t.profile) = &c.profile;
		}
		profile_time_ += timer.nanoseconds();
	}

	ExtensionPipeline::ExtensionQueue& queue_;
	const std::vector<ExtensionPipeline::Extension> batch_;
	vector<DP::AnchoredSwipe::Target<int16_t>> batch_targets_;
	size_t batch_pos_;
	int shard_;                            // shard of the queue popped from last
	bool queue_drained_ = false;           // the queue had no extension to pop, it is not asked again
	int64_t popped_;
	const DP::AnchoredSwipe::Config& cfg_;
	const ExtensionPipeline::ExtensionCallback& callback_;
	const int64_t padding_;
	alignas(32) int8_t score_table_[32 * 32];
	vector<Channel> channels_;
	vector<Done> done_;
	vector<ExtensionPipeline::Extension> out_;
	int64_t profile_time_, output_time_;
	const bool record_bands_;              // --band-hist is set
	vector<int64_t> band_hist_;            // bands of the extensions loaded, merged into the histogram of the run when done
};

/* Computes the extensions of the queue. Every extension is an independent local alignment of the
   query and the target of its pair, starting at the anchor and running into its own direction, and
   only its score and extent are reported. Merging the extensions of a pair into one alignment and
   filtering it is up to the callback.
   The batch, which the caller has popped from the queue, is computed first. Once it is used up, every
   channel of the kernel that completes pops the next extension from the queue, so extensions submitted
   while the kernel is running are computed by the same call. It returns once the queue had no extension
   for a channel and the channels loaded before are complete (see ExtensionSource). */
void anchored_swipe(std::vector<ExtensionPipeline::Extension>& batch, ExtensionPipeline::ExtensionQueue& queue, const DP::AnchoredSwipe::Config& cfg,
	const ExtensionPipeline::ExtensionCallback& callback)
{
	MEM_SCOPE("dp/anchored-swipe");
	TaskTimer total;

	if (batch.empty())
		return;

	{
		ExtensionSource source(queue, std::move(batch), cfg, callback);
		const DP::AnchoredSwipe::Options options = source.options();

		TaskTimer timer;
		const DP::AnchoredSwipe::Stats stats = DP::AnchoredSwipe::DISPATCH_ARCH::smith_waterman<::DISPATCH_ARCH::ScoreVector<int16_t, 0>>(source, cfg.band_max, options);
#ifdef ANCHORED_SWIPE_STAT
		cfg.stats.inc(Statistics::ANCHORED_LANE_CELLS, stats.lane_cells);
		cfg.stats.inc(Statistics::ANCHORED_LANE_CELLS_IDLE, stats.idle_cells);
		cfg.stats.inc(Statistics::ANCHORED_LANE_CELLS_TARGET_OOB, stats.target_oob_cells);
		cfg.stats.inc(Statistics::ANCHORED_LANE_CELLS_QUERY_OOB, stats.query_oob_cells);
		cfg.stats.inc(Statistics::ANCHORED_LANE_CELLS_OUT_OF_BAND, stats.out_of_band_cells);
#endif
		cfg.stats.inc(Statistics::SWIPE_TASKS_TOTAL);
		cfg.stats.inc(Statistics::TIME_SW, timer.microseconds());

		source.flush();
	}
	cfg.stats.inc(Statistics::TIME_ANCHORED_SWIPE, total.microseconds());
}

}

DISPATCH_4(void, anchored_swipe, std::vector<ExtensionPipeline::Extension>&, batch, ExtensionPipeline::ExtensionQueue&, queue, const DP::AnchoredSwipe::Config&, cfg,
	const ExtensionPipeline::ExtensionCallback&, callback)

}}
