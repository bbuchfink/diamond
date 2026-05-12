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

#include <set>
#include <algorithm>
#include <iomanip>
#include "taxonomy_nodes.h"
#include "util/log_stream.h"
#include "util/string/string.h"
#include "legacy/dmnd/io.h"
#include "blastdb/taxdmp.h"

using std::string;
using std::reverse;
using std::vector;
using std::set;
using std::endl;
using std::to_string;
using std::runtime_error;

TaxonomyNodes::TaxonomyNodes(const string& file_name)
{
	auto f = [&](TaxId taxid, TaxId parent, const string& rank) {
		parent_.resize(taxid + 1, 0);
		parent_[taxid] = parent;
		rank_.resize(taxid + 1, Rank::none);
		rank_[taxid] = Rank(rank.c_str());
		};
	read_nodes_dmp(file_name, f);
}

void TaxonomyNodes::save(Serializer &out)
{
	TaskTimer timer("Building taxonomy nodes");
	serialize(out, parent_);
	out.write(rank_.data(), rank_.size());
	timer.finish();
	message_stream << parent_.size() << " taxonomy nodes processed." << endl;
	size_t rank_count[Rank::count];
	std::fill(rank_count, rank_count + Rank::count, 0);
	for (const Rank r : rank_) {
		++rank_count[r];
	}
	
	const size_t w = MAX_LEN(Rank::names) + 2;
	message_stream << "Number of nodes assigned to rank:" << endl;
	for (size_t i = 0; i < Rank::count; ++i)
		message_stream << std::left << std::setw(w) << Rank::names[i] << rank_count[i] << endl;
	message_stream << endl;
}

TaxonomyNodes::TaxonomyNodes(Deserializer &in, uint32_t db_build)
{
	deserialize(in, parent_);
	if (db_build >= 131) {
		rank_.resize(parent_.size());
		in.read(rank_.data(), rank_.size());
	}
}