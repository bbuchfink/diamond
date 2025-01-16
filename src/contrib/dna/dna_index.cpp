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

#include <functional>
#include "dna_index.h"
#include "../data/sequence_set.h"
#include "../data/queries.h"
#include <atomic>
#include "../util/algo/partition.h"
#include <thread>
#include <queue>
#include <numeric>
#include "../util/util.h"
//#include "google/protobuf/arena.h"

using std::atomic;
using std::thread;
using std::vector;
using std::atomic_size_t;
using std::pair;



namespace Dna {

    Index::Index(Search::Config &cfg, char *ref_buffer) :
        ref_buffer_(ref_buffer) {
    const SeedPartitionRange range(0, Const::seedp);

    const SeedHistogram &ref_hst = cfg.target->hst();
    TaskTimer timer("Building reference seed array", true);


    const EnumCfg enum_ref{&ref_hst.partition(), 0, 1, cfg.seed_encoding, nullptr, false, false,
                           cfg.seed_complexity_cut,
                           MaskingAlgo::NONE, cfg.minimizer_window, false, false, cfg.sketch_size };


    seed_arr_.reset(new SeedArray(*cfg.target, ref_hst.get(0), range, ref_buffer, &no_filter, enum_ref));
    timer.go("Building reference index");
    count_minimizers(range);
    build_index(range,filter_repetitive(range));
}

Index::~Index() { delete[] ref_buffer_; }

pair<SeedArray::Entry *, SeedArray::Entry *> Index::contains(PackedSeed seed) const {
    unsigned partition = seed_partition(seed);
    unsigned key = seed_partition_offset(seed);

    if (dna_index_[partition]->size() == 0) {
        return {nullptr, nullptr};
    }

    auto hash_lookup = dna_index_[partition]->find_entry(key);

    if (hash_lookup == nullptr) {
        return {nullptr, nullptr};
    }
    SeedArray::Entry *first = seed_arr_->begin(partition) + hash_lookup->value; //
    SeedArray::Entry *end = seed_arr_->begin(partition) + seed_arr_->size(partition);
    for (auto i = first + 1; i != end; ++i) {
        if (i->key != first->key)
            return {first, i};
    }
    return {first, end};
}
void Index::count_worker(std::atomic<unsigned int> *seedp) {
    unsigned part;
    while ((part = (*seedp)++) < Const::seedp) {
        std::sort(seed_arr_->begin(part), seed_arr_->begin(part) + seed_arr_->size(part));
        auto it = merge_keys(seed_arr_->begin(part), (seed_arr_->begin(part) + seed_arr_->size(part)),
                             SeedArray::Entry::GetKey());

        unsigned count = 0;
        while (it.good()) {
            ++count;
            ++it;
        }
        this->minimizer_counts[part] = count;
    }
}
void Index::count_minimizers(const SeedPartitionRange &range) {
    std::function<void(atomic<unsigned> *seedp)> functor = [this](atomic<unsigned int> *PH1) { count_worker(PH1); };
    atomic<unsigned> seedp(range.begin());
    vector<std::thread> threads;
    for (int i = 0; i < config.threads_; ++i)
        threads.emplace_back(functor, &seedp);
    for (auto &t: threads)
        t.join();
    this->n_minimizer = std::accumulate(minimizer_counts.begin(), minimizer_counts.end(),static_cast<unsigned>(0));
}

void Index::index_worker(atomic<unsigned> *seedp, int cutoff) {
    unsigned part;
    while ((part = (*seedp)++) < Const::seedp) {
        dna_index_[part].reset(new HashTable<SeedOffset, unsigned, MurmurHash, Modulo>(this->minimizer_counts[part] * 1.2, MurmurHash()));

        auto it2 = merge_keys(seed_arr_->begin(part), (seed_arr_->begin(part) + seed_arr_->size(part)),SeedArray::Entry::GetKey());
        while (it2.good()) {
           if (it2.count() < cutoff)
                dna_index_[part]->insert(it2.key())->value = it2.begin() - seed_arr_->begin(part);
            ++it2;
        }
    }
}


void Index::build_index(const SeedPartitionRange &range, int repetitive_cutoff) {

    std::function<void(atomic<unsigned> *seedp, int cutoff)> functor = [this](atomic<unsigned int> *PH1, int cutoff) { index_worker(PH1,cutoff); };
    atomic<unsigned> seedp(range.begin());
    vector<std::thread> threads;
    for (int i = 0; i < config.threads_; ++i)
        threads.emplace_back(functor, &seedp, repetitive_cutoff);
    for (auto &t: threads)
        t.join();
}


void Index::filter_worker(std::atomic<unsigned> *seedp, int index, std::vector<std::priority_queue<int, std::vector<int>, std::greater<>>> &rep_thread, int n) {

    unsigned part;
    while ((part = (*seedp)++) < Const::seedp) {
        auto it = merge_keys(seed_arr_->begin(part), (seed_arr_->begin(part) + seed_arr_->size(part)),
                             SeedArray::Entry::GetKey());
        while (it.good()) {
            if (rep_thread[index].size() < n) {
                rep_thread[index].emplace(it.count());
            } else if (rep_thread[index].top() < it.count()){
                rep_thread[index].emplace(it.count());
                rep_thread[index].pop();
            }
            ++it;
        }
    }
}

int Index::filter_repetitive(const SeedPartitionRange &range) {

    const int n = n_minimizer * config.repetitive_cutoff;

    if(n < 1)
        return INT_MAX;


    std::vector<std::priority_queue<int, std::vector<int>, std::greater<>>> rep_thread(config.threads_);

    std::function<void(atomic<unsigned> *seedp, int heap_index)> functor = [this, &rep_thread, n](atomic<unsigned int> *PH1, int heap_index) { filter_worker(PH1,heap_index, rep_thread, n); };
    std::priority_queue<int,std::vector<int>,std::greater<>> repetitive_max;

    atomic<unsigned> seedp(range.begin());
    vector<std::thread> threads;

    for (int i = 0; i < config.threads_; ++i) {
        threads.emplace_back(functor, &seedp, i);
    }
    for (auto &t: threads)
        t.join();

    for(auto &heap: rep_thread){
        while (!heap.empty()) {
            if (repetitive_max.size() < n) {
                repetitive_max.emplace(heap.top());
            } else if (heap.top() > repetitive_max.top()) {
                repetitive_max.emplace(heap.top());
                repetitive_max.pop();
            }
            heap.pop();
        }
    }
    return repetitive_max.top();
}

}