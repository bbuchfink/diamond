#include <unordered_map>
#include <array>
#include <iostream>
#include <vector>
#include <fstream>
#include "tsv_record.h"
#include "../basic/config.h"
#include "../basic/value.h"
#include "../util/util.h"
#include "../basic/sequence.h"
#include "../util/seq_file_format.h"
#include "../util/sequence/sequence.h"
#include "../stats/cbs.h"
#include "../util/algo/MurmurHash3.h"
#include "../util/sequence/sequence.h"
#include "../data/sequence_file.h"
#include "../search/search.h"
#include "../data/enum_seeds.h"
#define _REENTRANT
#include "../lib/ips4o/ips4o.hpp"
#include "../masking/masking.h"
#include "../util/ptr_vector.h"
#include "../util/string/tokenizer.h"
#include "../util/algo/algo.h"
#include "../util/data_structures/hash_table.h"
#include "../util/string/fixed_string.h"
#include "../util/tsv/table.h"
#include "../util/tsv/file.h"

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
using Util::Tsv::TokenIterator;
using Util::Seq::FastaIterator;

void filter_blasttab() {
	TextInputFile in("");
	TSVRecord r;
	string query;
	size_t query_hit;
	while (in >> r) {
		if (r.qseqid != query) {
			query = r.qseqid;
			query_hit = 0;
		}
		else
			++query_hit;
		/*if (query_hit < config.max_alignments && r.evalue <= config.max_evalue)
			cout << r << endl;*/
	}
}

void split() {
	TextInputFile in(config.single_query_file());
	string id;
	vector<Letter> seq;
	size_t n = 0, f = 0, b = (size_t)(config.chunk_size * 1e9), seqs = 0;
	OutputFile *out = new OutputFile(std::to_string(f) + ".faa.zst", Compressor::ZSTD);
	while (FASTA_format().get_seq(id, seq, in, value_traits)) {
		if (n >= b) {
			out->close();
			delete out;
			out = new OutputFile(std::to_string(++f) + ".faa.zst", Compressor::ZSTD);
			n = 0;
		}
		string blast_id = Util::Seq::seqid(id.c_str(), false);
		Util::Seq::format(Sequence(seq), blast_id.c_str(), nullptr, *out, "fasta", amino_acid_traits);
		n += seq.size();
		++seqs;
		if (seqs % 1000000 == 0)
			std::cerr << "#Sequences processed: " << seqs << " #letters:" << n << endl;
	}
	out->close();
	delete out;
}

void composition() {
	TextInputFile in(config.single_query_file());
	string id;
	vector<Letter> seq;
	while (FASTA_format().get_seq(id, seq, in, value_traits)) {
		auto c = Stats::composition(seq);
		for (double x : c)
			std::cout << x << '\t';
		std::cout << endl;
	}
}

void hash_seqs() {
	TextInputFile f(config.query_file.front());
	FASTA_format fmt;
	string id;
	vector<Letter> seq;
	while (fmt.get_seq(id, seq, f, amino_acid_traits)) {
		array<char, 16> hash;
		hash.fill('\0');
		MurmurHash3_x64_128(seq.data(), (int)seq.size(), hash.data(), hash.data());
		cout << Util::Seq::seqid(id.c_str(), false) << '\t' << hex_print(hash.data(), 16) << endl;
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
	::shapes = ShapeConfig(config.shape_mask.empty() ? Search::shape_codes[(int)align_mode.sequence_type].at(Sensitivity::DEFAULT) : config.shape_mask, config.shapes);
	Reduction::reduction = Reduction("A R N D C Q E G H I L K M F P S T W Y V");
	EnumCfg cfg{ &parts, 0, 1, SeedEncoding::SPACED_FACTOR, nullptr, false, false, config.seed_cut_, MaskingAlgo::NONE, 0 };
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
		s = Reduction::reduction.decode_seed(i->second, shapes[0].weight_);
		cout << i->first << '\t' << s << endl;
	}
}

using Acc = FixedString<32>;

void index_fasta() {
	HashTable<Acc, int64_t, Acc::Hash, Modulo> table(size_t(config.max_target_seqs_ * 1.2), Acc::Hash());
	std::ifstream f(config.query_file.front());
	string l;
	int64_t n = 0;
	for (;;) {
		int64_t p = f.tellg();
		if (!std::getline(f, l))
			break;
		if (l.empty() || l[0] != '>')
			continue;
		auto e = table.insert(Acc(Util::Seq::seqid(l.c_str() + 1, false)));
		e->value = p + 1;
		++n;
	}
	OutputFile out(config.query_file.front() + ".htidx");
	out.write(table.data(), table.size());
	out.close();
	message_stream << "#Sequences: " << n << endl;
}

void fetch_seq() {

}

void length_sort() {
	Schema schema{ Type::STRING, Type::STRING };
	Util::Tsv::Table table(schema);
}

#ifdef EXTRA

void sort() {
	Util::Tsv::File input(Schema{ Type::INT64, Type::STRING }, config.query_file.front());
	//auto f = input.sort(0, 1, config.output_file);
	//delete f;
}

#endif