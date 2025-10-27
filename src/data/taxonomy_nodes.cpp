/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

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