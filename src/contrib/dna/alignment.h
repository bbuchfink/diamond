/****
DIAMOND protein aligner
Copyright (C) 2023 Vincent Spath

Code developed by Vincent Spath <vincent.spath@tuebingen.mpg.de>

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

#include "../data/block/block.h"
#include "../align/extend.h"
#include "../lib/ksw2/ksw2.h"
#include "seed_set_dna.h"
#include "../util/util.h"
#include "../util/sequence/translate.h"
#include "../basic/shape_config.h"
#include <iostream>


#pragma once


namespace Dna{
const int KSW2_ZDROP_EXTENSION = 40; // TODO: should be different for short and long reads
const int KSW2_ZDROP_BETWEEN_ANCHORS = 100;
const int KSW2_BAND_EXTENSION = 40;
const int KSW2_BAND_GLOBAL = 30;
const int WFA_BAND_EXTENSION = 20;
const int WFA_ZDROP_EXTENSION = 100;
const int WFA_ZDROP_GLOBAL = 500;

    struct Cigar {

    public:
        Cigar() : query_extension_distance_(0), target_extension_distance_(0) {}

        explicit Cigar(std::size_t reserveSize = 100) : query_extension_distance_(0), target_extension_distance_(0) {
            cigar_data.reserve(reserveSize);
        }

        void reserveCigarSpace(std::size_t reserveSize) {
            cigar_data.reserve(reserveSize);
        }

        void extendCigar(const std::vector<std::pair<int, char>>& otherVector) {
            this->cigar_data.insert(this->cigar_data.end(), otherVector.begin(), otherVector.end());
        }

        void extendCigar(uint32_t length, char cigarOperation) {
            cigar_data.emplace_back(length, cigarOperation);
        }

        int query_extension_distance() const { return query_extension_distance_; }

        int target_extension_distance() const { return target_extension_distance_; }


        void setMaxValues(int queryStart, int targetStart) {
            query_extension_distance_ = queryStart;
            target_extension_distance_ = targetStart;
        }

        const std::vector<std::pair<int, char>>& getCigarDataConst() const {
            return cigar_data;
        }

        std::vector<std::pair<int, char>>& getCigarData() {
            return cigar_data;
        }

        int score = 0;
        int peakScore = 0;
        int peakScoreCigarIndex = 0;
        int peakScoreAnchorIndex = 0;

    private:
        int query_extension_distance_{0}, target_extension_distance_{0};
        std::vector<std::pair<int, char>> cigar_data;
    };

    enum AlignmentStatus : uint8_t {
        NOT_DROPPED = 0,
        DROPPED = 1,
        NEGATIVE_SCORE = 2
    };



    AlignmentStatus compute_ksw_cigar(const Sequence &targetSequence, const Sequence &querySequence, const Search::Config &cfg,
                                      int flag,
                                      Cigar &extension, const int zdrop, const int band);
    AlignmentStatus
    compute_wfa_cigar(const Search::Config &cfg, const std::string &querySequence, Cigar &extension, bool left,
                      bool global,
                      const std::string &targetSequence, const int band);

    Hsp build_hsp_from_cigar(const Cigar &cigar, const Sequence &target, const Sequence &query, const int firstAnchor_i, const int firstAnchor_j, bool isReverse, const Search::Config &cfg);
}
