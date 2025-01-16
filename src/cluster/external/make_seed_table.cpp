/****
DIAMOND protein aligner
Copyright (C) 2024 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <buchfink@gmail.com>

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

#include "basic/config.h"
#include "data/sequence_file.h"
#include "data/enum_seeds.h"
#include "search/search.h"
#include "util/tsv/tsv.h"
#include "util/log_stream.h"

using std::thread;
using std::unique_ptr;
using std::endl;
using std::vector;
using namespace Util::Tsv;

namespace Cluster {

struct SeedCallback {
	SeedCallback(File& out, const Block& block):
		out(out),
		block(block)
	{}
	bool operator()(uint64_t seed, uint64_t pos, uint32_t block_id, size_t shape) const {
		out.write_record(seed, block.block_id2oid(block_id), block.seqs().length(block_id));
		return true;
	};
	void finish() {
	}
	File& out;
	const Block& block;
};

void make_seed_table() {
	static const int64_t BLOCK_SIZE = 1000000000;
	config.database.require();
	const Sensitivity sens = Sensitivity::FASTER;
	
	::shapes = ShapeConfig(config.shape_mask.empty() ? Search::shape_codes.at(sens) : config.shape_mask, config.shapes);

	File out(Schema{ Type::INT64, Type::INT64, Type::INT64 }, config.output_file, Flags::WRITE);

	auto worker = [&out, sens] {
		TaskTimer timer("Opening the database");
		unique_ptr<SequenceFile> db(SequenceFile::auto_create({ config.database }, SequenceFile::Flags(), SequenceFile::Metadata(), value_traits));
		timer.finish();
		message_stream << "Sequences: " << db->sequence_count() << endl;
		message_stream << "Letters: " << db->letters() << endl;
		unique_ptr<Block> block;
		int n = 0;
		for (;;) {
			timer.go("Loading sequences");
			block.reset(db->load_seqs(BLOCK_SIZE));
			if (block->empty())
				break;
			timer.finish();
			message_stream << "Processing block " << (n++) + 1 << endl;
			const vector<BlockId> partition = block->seqs().partition(1);
			PtrVector<SeedCallback> cb;
			cb.push_back(new SeedCallback(out, *block));
			EnumCfg cfg{
				&partition,
				0,
				1,
				SeedEncoding::SPACED_FACTOR,
				nullptr,
				false,
				false,
				0,
				MaskingAlgo::TANTAN | MaskingAlgo::MOTIF,
				Search::sensitivity_traits.at(sens).minimizer_window,
				false,
				false,
				Search::sensitivity_traits.at(sens).sketch_size
			};
			enum_seeds<SeedCallback, NoFilter>(*block, cb, nullptr, cfg);
		}
	};

	thread t(worker);
	t.join();
}

}