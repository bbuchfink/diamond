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

#include "../lib/ksw2/ksw2.h"
#include "../data/block/block.h"
#include "../align/extend.h"

#pragma once

const int KSW2_END_BONUS = 100;

namespace Dna {
/**
struct ExtensionTimer {
public:
    ExtensionTimer();
    void update(int operation , std::chrono::duration<long, std::ratio<1, 1000000000>> duration);
    std::chrono::time_point<std::chrono::system_clock, std::chrono::duration<long, std::ratio<1, 1000000000>>> total_time;
    std::chrono::time_point<std::chrono::system_clock, std::chrono::duration<long, std::ratio<1, 1000000000>>> seed_lookup_time;
    std::chrono::time_point<std::chrono::system_clock, std::chrono::duration<long, std::ratio<1, 1000000000>>> ungapped;
    std::chrono::time_point<std::chrono::system_clock, std::chrono::duration<long, std::ratio<1, 1000000000>>> trimming;
    std::chrono::time_point<std::chrono::system_clock, std::chrono::duration<long, std::ratio<1, 1000000000>>> ksw_extension;
    std::chrono::time_point<std::chrono::system_clock, std::chrono::duration<long, std::ratio<1, 1000000000>>> diamond_cigar;
};
 **/

    std::pair<std::vector<Extension::Match>, Extension::Stats> extend(const Search::Config &cfg, const Sequence &query);//, ExtensionTimer& extend_time);
}



