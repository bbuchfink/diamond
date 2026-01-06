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

#include <stdexcept>
#include "basic/config.h"
#include "search/search.h"
#include "seed_set.h"
#include "dmnd/dmnd.h"
#include "basic/shape_config.h"
#include "util/log_stream.h"

void makeindex() {
	static const size_t MAX_LETTERS = 100000000;
	if (config.database.empty())
		throw std::runtime_error("Missing parameter: database file (--db/-d).");
	DatabaseFile db(config.database);
	if (db.ref_header.letters > MAX_LETTERS)
		throw std::runtime_error("Indexing is only supported for databases of < 100000000 letters.");

	::shapes = ShapeConfig(config.shape_mask.empty() ? Search::shape_codes.at(config.sensitivity) : config.shape_mask, config.shapes);
	config.algo = Config::Algo::DOUBLE_INDEXED;

	db.flags() |= SequenceFile::Flags::SEQS;
	Block* block = db.load_seqs(MAX_LETTERS, nullptr);

	TaskTimer timer("Building index");
	HashedSeedSet index(*block, nullptr, 0.0, Search::soft_masking_algo(Search::sensitivity_traits.at(config.sensitivity)));

	timer.go("Writing to disk");
	OutputFile out(db.file_name() + ".seed_idx");
	out.write(SEED_INDEX_MAGIC_NUMBER);
	out.write(SEED_INDEX_VERSION);
	out.write((uint32_t)shapes.count());

	for (int i = 0; i < shapes.count(); ++i)
		out.write(index.table(i).size());

	for (int i = 0; i < shapes.count(); ++i) {
		out.write(index.table(i).data(), index.table(i).size() + HashedSeedSet::Table::PADDING);
	}

	out.close();
	db.close();
	delete block;
}