#include <stdexcept>
#include "basic/config.h"
#include "search/search.h"
#include "seed_set.h"
#include "dmnd/dmnd.h"
#include "basic/shape_config.h"
#include "util/log_stream.h"

void makeindex() {
	static const size_t MAX_LETTERS = 100'000'000;
	if (config.database.empty())
		throw std::runtime_error("Missing parameter: database file (--db/-d).");
	DatabaseFile db(config.database);
	if (db.ref_header.letters > MAX_LETTERS)
		throw std::runtime_error("Indexing is only supported for databases of < 100000000 letters.");

	::shapes = ShapeConfig(config.shape_mask.empty() ? Search::shape_codes.at(config.sensitivity) : config.shape_mask, config.shapes);
	config.algo = Config::Algo::DOUBLE_INDEXED;

	Block* block = db.load_seqs(MAX_LETTERS, nullptr, SequenceFile::LoadFlags::SEQS);

	TaskTimer timer("Building index");
	HashedSeedSet index(*block, nullptr, 0.0, Search::soft_masking_algo(Search::sensitivity_traits.at(config.sensitivity)));

	timer.go("Writing to disk");
	OutputFile out(db.file_name() + ".seed_idx");
	out.write(SEED_INDEX_MAGIC_NUMBER);
	out.write(SEED_INDEX_VERSION);
	out.write((uint32_t)shapes.count());

	for (unsigned i = 0; i < shapes.count(); ++i)
		out.write(index.table(i).size());

	for (unsigned i = 0; i < shapes.count(); ++i) {
		out.write(index.table(i).data(), index.table(i).size() + HashedSeedSet::Table::PADDING);
	}

	out.close();
	db.close();
	delete block;
}