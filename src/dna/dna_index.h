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

#include <atomic>
#include "../data/seed_array.h"
#include "../util/data_structures/hash_table.h"

#pragma once


namespace Dna{

class Index{

public:
    Index(Search::Config& cfg,char *ref_buffer);
    std::pair<SeedArray::Entry*,SeedArray::Entry*> contains(PackedSeed seed) const;
    ~Index();




private:
    void index_worker(std::atomic<unsigned> *seedp);
    void build_index(const SeedPartitionRange& range);
    std::unique_ptr<SeedArray> seed_arr_;
    std::array<std::unique_ptr<HashTable<SeedOffset,unsigned,MurmurHash,Modulo>>, Const::seedp> dna_index_;
    char* ref_buffer_;


};

}