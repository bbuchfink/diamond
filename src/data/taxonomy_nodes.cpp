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

const std::map<char, char> Rank::legacy_rank_codes = {
    {0, Rank::none},
    {1, Rank::superkingdom},
    {2, Rank::kingdom},
    {3, Rank::subkingdom},
    {4, Rank::superphylum},
    {5, Rank::phylum},
    {6, Rank::subphylum},
    {7, Rank::superclass},
    {8, Rank::class_rank},
    {9, Rank::subclass},
    {10, Rank::infraclass},
    {11, Rank::cohort},
    {12, Rank::subcohort},
    {13, Rank::superorder},
    {14, Rank::order},
    {15, Rank::suborder},
    {16, Rank::infraorder},
    {17, Rank::parvorder},
    {18, Rank::superfamily},
    {19, Rank::family},
    {20, Rank::subfamily},
    {21, Rank::tribe},
    {22, Rank::subtribe},
    {23, Rank::genus},
    {24, Rank::subgenus},
    {25, Rank::section},
    {26, Rank::subsection},
    {27, Rank::series},
    {28, Rank::species_group},
    {29, Rank::species_subgroup},
    {30, Rank::species},
    {31, Rank::subspecies},
    {32, Rank::varietas},
    {33, Rank::forma},
    {34, Rank::strain},
    {35, Rank::biotype},
    {36, Rank::clade},
    {37, Rank::forma_specialis},
    {38, Rank::genotype},
    {39, Rank::isolate},
    {40, Rank::morph},
    {41, Rank::pathogroup},
    {42, Rank::serogroup},
    {43, Rank::serotype},
    {44, Rank::subvariety}
};

TaxonomyNodes::TaxonomyNodes(const string& file_name):
	node_count_(0)
{
	auto f = [&](TaxId taxid, TaxId parent, const string& rank) {
		if (safe_cast<size_t>(taxid) >= parent_.size())
			parent_.resize(taxid + 1, 0);
		parent_[taxid] = parent;
		if (safe_cast<size_t>(taxid) >= rank_.size())
			rank_.resize(taxid + 1, Rank::none);
		rank_[taxid] = Rank(rank.c_str());
		++node_count_;
		};
	read_nodes_dmp(file_name, f);
}

void TaxonomyNodes::save(Serializer &out)
{
	TaskTimer timer("Building taxonomy nodes");
	serialize(out, parent_);
	out.write(rank_.data(), rank_.size());
	timer.finish();
    *message_stream << node_count_ << " taxonomy nodes processed." << endl;
    *message_stream << "Maximum taxon id: " << parent_.size() - 1 << endl;
	size_t rank_count[Rank::count];
	std::fill(rank_count, rank_count + Rank::count, 0);
	for (const Rank r : rank_) {
		++rank_count[r];
	}
	
	const size_t w = MAX_LEN(Rank::names) + 2;
    *message_stream << "Number of nodes assigned to rank:" << endl;
	for (size_t i = 0; i < Rank::count; ++i)
        *message_stream << std::left << std::setw(w) << Rank::names[i] << rank_count[i] << endl;
    *message_stream << endl;
}

TaxonomyNodes::TaxonomyNodes(File &in, uint32_t db_build):
	node_count_(0)
{
	deserialize(in, parent_);
	if (db_build >= 131) {
		rank_.resize(parent_.size());
		in.read(rank_.data(), rank_.size());
	}
}