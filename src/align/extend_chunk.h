static size_t lazy_masking(vector<uint32_t>::const_iterator target_block_ids, vector<uint32_t>::const_iterator target_block_ids_end, Block& targets, const MaskingAlgo algo) {
	if (algo == MaskingAlgo::NONE)
		return 0;
	vector<Letter> seq;
	const Masking& masking = Masking::get();
	size_t n = 0;
	for (auto i = target_block_ids; i != target_block_ids_end; ++i)
		if (targets.fetch_seq_if_unmasked(*i, seq)) {
			masking(seq.data(), seq.size(), algo, *i);
			targets.write_masked_seq(*i, seq);
			++n;
		}
	return n;
}

pair<vector<Target>, Stats> extend(BlockId query_id,
	const Sequence *query_seq,
	int source_query_len,
	const Bias_correction *query_cb,
	const ::Stats::Composition& query_comp,
	FlatArray<SeedHit>::Iterator seed_hits,
	FlatArray<SeedHit>::Iterator seed_hits_end,
	vector<uint32_t>::const_iterator target_block_ids,
	const Search::Config& cfg,
	Statistics& stat,
	DP::Flags flags,
	const HspValues hsp_values)
{
	static const Loc GAPPED_FILTER_MIN_QLEN = 85;
	const int64_t n = seed_hits_end - seed_hits;
	stat.inc(Statistics::TARGET_HITS2, n);
	task_timer timer(flag_any(flags, DP::Flags::PARALLEL) ? config.target_parallel_verbosity : UINT_MAX);

	if (cfg.lazy_masking && !config.global_ranking_targets)
		stat.inc(Statistics::MASKED_LAZY, lazy_masking(target_block_ids, target_block_ids + n, *cfg.target, cfg.target_masking));

	pair<FlatArray<SeedHit>, vector<uint32_t>> gf;
	if (cfg.gapped_filter_evalue > 0.0 && config.global_ranking_targets == 0 && (!align_mode.query_translated || query_seq[0].length() >= GAPPED_FILTER_MIN_QLEN)) {
		timer.go("Computing gapped filter");
		gf = gapped_filter(query_seq, query_cb, seed_hits, seed_hits_end, target_block_ids, stat, flags, cfg);
		if (!flag_any(flags, DP::Flags::PARALLEL))
			stat.inc(Statistics::TIME_GAPPED_FILTER, timer.microseconds());
		seed_hits = gf.first.begin();
		seed_hits_end = gf.first.end();
		target_block_ids = gf.second.cbegin();
	}
	stat.inc(Statistics::TARGET_HITS3, seed_hits_end - seed_hits);

	timer.go("Computing chaining");
	vector<WorkTarget> targets = ungapped_stage(query_seq, query_cb, query_comp, seed_hits, seed_hits_end, target_block_ids, flags, stat, *cfg.target, cfg.extension_mode);
	if (!flag_any(flags, DP::Flags::PARALLEL))
		stat.inc(Statistics::TIME_CHAINING, timer.microseconds());

	return align(targets, query_seq, cfg.query->ids()[query_id], query_cb, source_query_len, flags, hsp_values, cfg.extension_mode, *cfg.thread_pool, cfg, stat);
}
