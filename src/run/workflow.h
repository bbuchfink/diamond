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

#pragma once
#include <memory>
#include <vector>
#include "data/sequence_file.h"
#include "util/io/file.h"

void get_seq();
void random_seqs();

namespace Search {

void run(std::unique_ptr<std::vector<BitVector>>& target_seed_hits, const std::shared_ptr<SequenceFile>& db = nullptr, const std::shared_ptr<SequenceFile>& query = nullptr, const std::shared_ptr<File>& out = nullptr, const std::shared_ptr<DbFilter>& db_filter = nullptr);

}
