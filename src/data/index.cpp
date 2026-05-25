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
	Block* block = db.load_seqs(MAX_LETTERS, 0, nullptr);

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