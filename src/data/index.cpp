#include <stdexcept>
#include "../basic/config.h"
#include "reference.h"
#include "seed_set.h"
#include "../search/search.h"
#include "seed_set.h"
#include "dmnd/dmnd.h"

void makeindex() {
	static const size_t MAX_LETTERS = 100000000;
	if (config.database.empty())
		throw std::runtime_error("Missing parameter: database file (--db/-d).");
	DatabaseFile db(config.database);
	if (db.ref_header.letters > MAX_LETTERS)
		throw std::runtime_error("Indexing is only supported for databases of < 100000000 letters.");

	setup_search();
	config.algo = Config::Algo::DOUBLE_INDEXED;

	vector<uint32_t> block2db_id;
	SequenceSet* seqs;
	db.load_seqs(&block2db_id, MAX_LETTERS, &seqs, nullptr, false);

	task_timer timer("Building index");
	HashedSeedSet index(*seqs);

	timer.go("Writing to disk");
	OutputFile out(db.file_name() + ".seed_idx");
	out.write(SEED_INDEX_MAGIC_NUMBER);
	out.write(SEED_INDEX_VERSION);
	out.write((uint32_t)shapes.count());

	for (unsigned i = 0; i < shapes.count(); ++i)
		out.write(index.table(i).size());

	for (unsigned i = 0; i < shapes.count(); ++i) {
		out.write(index.table(i).data(), index.table(i).size() + HashedSeedSet::HashSet::PADDING);
	}

	out.close();
	db.close();
	delete seqs;
}