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

#include <unordered_map>
#include <array>
#include <iostream>
#include <vector>
#include "basic/config.h"
#include "basic/value.h"
#include "util/util.h"
#include "basic/sequence.h"
#include "util/sequence/sequence.h"
#include "murmurhash/MurmurHash3.h"
#include "data/sequence_file.h"
#include "search/search.h"
#include "data/enum_seeds.h"
#define _REENTRANT
#include "ips4o/ips4o.hpp"
#include "masking/masking.h"
#include "util/ptr_vector.h"
#include "util/tsv/table.h"
#include "util/tsv/tsv.h"
#include "util/string/tokenizer_dyn.h"

using std::array;
using std::cout;
using std::endl;
using std::vector;
using std::unique_ptr;
using std::pair;
using std::unordered_map;
using std::cerr;
using std::ofstream;
using std::string;
using Util::Tsv::Schema;
using Util::Tsv::Type;
using Util::String::TokenIterator;
using std::runtime_error;

void split() {
	unique_ptr<SequenceFile> in;
	try {
		in.reset(SequenceFile::auto_create({ config.single_query_file() }));
	}
	catch (FormatDetectionError& e) {
		throw runtime_error("Error opening input file: " + std::string(e.what()));
	}
	string id;
	vector<Letter> seq;
	size_t n = 0, f = 0, b = (size_t)(config.chunk_size * 1e9), seqs = 0;
	OutputFile *out = new OutputFile(std::to_string(f) + ".faa.zst", Compressor::ZSTD);
	TextBuffer buf;
	while (in->read_seq(seq, id, nullptr)) {
		if (n >= b) {
			out->close();
			delete out;
			out = new OutputFile(std::to_string(++f) + ".faa.zst", Compressor::ZSTD);
			n = 0;
		}
		string blast_id = Util::Seq::seqid(id.c_str());
		Util::Seq::format(Sequence(seq), blast_id.c_str(), nullptr, buf, "fasta", amino_acid_traits);
		out->write(buf.data(), buf.size());
		buf.clear();
		n += seq.size();
		++seqs;
		if (seqs % 1000000 == 0)
			std::cerr << "#Sequences processed: " << seqs << " #letters:" << n << endl;
	}
	out->close();
	delete out;
}

void hash_seqs() {
	TextInputFile f(config.query_file.front());
	//FASTA_format fmt;
	string id;
	vector<Letter> seq;
	//while (fmt.get_seq(id, seq, f, amino_acid_traits)) {
	while(true) {
		array<char, 16> hash;
		hash.fill('\0');
		MurmurHash3_x64_128(seq.data(), (int)seq.size(), hash.data(), hash.data());
		cout << Util::Seq::seqid(id.c_str()) << '\t' << hex_print(hash.data(), 16) << endl;
	}
	f.close();
}

static double freq(const string& s, const Reduction& r) {
	double f = 0.0;
	for (char c : s) {
		f += r.freq(r(amino_acid_traits.from_char(c)));
	}
	return f / s.length();
}

void list_seeds() {
	struct Callback {
		bool operator()(uint64_t seed, size_t, unsigned, size_t) {
			seeds.push_back(seed);
			return true;
		};
		void finish() {}
		vector<uint64_t>& seeds;
	};
	unique_ptr<SequenceFile> db(SequenceFile::auto_create({ config.database }));
	unique_ptr<Block> block(db->load_seqs(INT64_MAX));
	mask_seqs(block->seqs(), Masking::get(), true, MaskingAlgo::TANTAN);
	vector<uint64_t> seeds;
	seeds.reserve(block->seqs().letters());
	PtrVector<Callback> cb;
	cb.push_back(new Callback{ seeds });
	auto parts = block->seqs().partition(1);
	::shapes = ShapeConfig(config.shape_mask.empty() ? Search::shape_codes.at(Sensitivity::DEFAULT) : config.shape_mask, config.shapes);
	Reduction custom("A R N D C Q E G H I L K M F P S T W Y V");
	Reduction::set_reduction(custom);
	EnumCfg cfg{ &parts, 0, 1, SeedEncoding::SPACED_FACTOR, nullptr, false, false, config.seed_cut_, MaskingAlgo::NONE, 0, false, false, 0, nullptr };
	enum_seeds(*block, cb, &no_filter, cfg);
	ips4o::parallel::sort(seeds.begin(), seeds.end());

	auto it = merge_keys(seeds.begin(), seeds.end(), [](uint64_t seed) {return seed; });
	vector<pair<uint64_t, uint64_t>> counts;
	while (it.good()) {
		counts.push_back({ it.count(), it.key() });
		++it;
	}
	ips4o::parallel::sort(counts.begin(), counts.end());

	auto end = std::min(counts.rbegin() + config.query_count, counts.rend());
	string s;
	for (auto i = counts.rbegin(); i != end; ++i) {
		s = Reduction::get_reduction().decode_seed(i->second, shapes[0].weight_);
		cout << i->first << '\t' << s << endl;
	}
}