/****
DIAMOND protein aligner
Copyright (C) 2022 Dimitrios Koutsogiannis

Code developed by Dimitrios Koutsogiannis <dimitrios.koutsogiannis@tue.mpg.de>

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
#include "dna_index.h"
#include "../data/sequence_set.h"
#include "../data/queries.h"
#include <atomic>
#include "../util/algo/partition.h"
#include <thread>
#include "../util/util.h"

using std::atomic;
using std::thread;
using std::vector;
using std::atomic_size_t;
using std::pair;


namespace Dna {

Index::Index(Search::Config &cfg, char *ref_buffer):
ref_buffer_(ref_buffer)
{
    const SeedPartitionRange range(0, Const::seedp);

    SequenceSet &ref_seqs = cfg.target->seqs();
    const SeedHistogram &ref_hst = cfg.target->hst();
    TaskTimer timer("Building reference seed array", true);


    const EnumCfg enum_ref{&ref_hst.partition(), 0, 1, cfg.seed_encoding, nullptr, false, false,
                           cfg.seed_complexity_cut,
                           MaskingAlgo::NONE, cfg.minimizer_window, false, false };


    seed_arr_.reset(new SeedArray(*cfg.target, ref_hst.get(0), range, ref_buffer, &no_filter, enum_ref));

    build_index(range);

}
Index::~Index(){delete [] ref_buffer_;}

pair<SeedArray::Entry*,SeedArray::Entry*>  Index::contains(PackedSeed seed) const {
    unsigned partition = seed_partition(seed);
    unsigned key = seed_partition_offset(seed);

    if(dna_index_[partition]->size() == 0){
        return {nullptr, nullptr};
    }

    auto hash_lookup = dna_index_[partition]->find_entry(key);

    if(hash_lookup == nullptr){
        return {nullptr, nullptr};
    }
    SeedArray::Entry* first =  seed_arr_->begin(partition) + hash_lookup->value; //
    SeedArray::Entry* end = seed_arr_->begin(partition) + seed_arr_->size(partition);
    for(auto i = first+1; i != end; ++i){
       if(i->key != first->key)
            return {first,i-1};
    }
    return {first,end};
}
void Index::index_worker(atomic<unsigned> *seedp) {
    unsigned part;
    while ((part = (*seedp)++) < Const::seedp) {
        std::sort(seed_arr_->begin(part), seed_arr_->begin(part) + seed_arr_->size(part));
        auto it = merge_keys(seed_arr_->begin(part), (seed_arr_->begin(part) + seed_arr_->size(part)), SeedArray::Entry::GetKey());
        unsigned count = 0;
        while(it.good()){
            ++count;
            ++it;
        }


        dna_index_[part].reset(new HashTable<SeedOffset, unsigned, MurmurHash, Modulo>(count * 1.2, MurmurHash()));
        auto it2 = merge_keys(seed_arr_->begin(part), (seed_arr_->begin(part) + seed_arr_->size(part)), SeedArray::Entry::GetKey());
        while(it2.good()){
           dna_index_[part]->insert(it2.key())->value = it2.begin() - seed_arr_->begin(part);
           ++it2;
        }
    }
}
void Index::build_index(const SeedPartitionRange& range){
    std::function<void (atomic <unsigned> *seedp)> functor = [this](atomic<unsigned int>* PH1) {index_worker(PH1);};
    atomic<unsigned> seedp(range.begin());
    vector<std::thread> threads;
    for (int i = 0; i < config.threads_; ++i)
        threads.emplace_back(functor, &seedp);
    for (auto &t: threads)
        t.join();
}
}
