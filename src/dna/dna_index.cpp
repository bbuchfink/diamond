#include "dna_index.h"
#include "../data/sequence_set.h"
#include "../data/queries.h"



namespace Dna{
Index::Index(Search::Config& cfg,char *ref_buffer)
{
        SequenceSet& ref_seqs = cfg.target->seqs();
        const SeedHistogram& ref_hst = cfg.target->hst();
        task_timer timer("Building reference seed array", true);


        const EnumCfg enum_ref{&ref_hst.partition(), 0, 1, cfg.seed_encoding, nullptr, false, false, cfg.seed_complexity_cut,
                                MaskingAlgo::NONE,cfg.minimizer_window };

        ref_idx_ =  std::make_unique<SeedArray>(*cfg.target, ref_hst.get(0), SeedPartitionRange(), ref_buffer, &no_filter, enum_ref);

}
}
