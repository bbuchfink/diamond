#include <mutex>
#include <unordered_map>
#include <thread>
#include "../basic/config.h"
#include "../data/sequence_file.h"
#include "../util/string/tokenizer.h"
#include "../util/tsv/tsv.h"
#include "../dp/dp.h"
#include "../output/output_format.h"
#include "../output/output.h"
#include "../masking/masking.h"
#include "../util/sequence/sequence.h"
#include "../util/string/fixed_string.h"

using std::unique_ptr;
using std::endl;
using std::mutex;
using std::lock_guard;
using std::thread;
using std::list;
using std::array;
using std::thread;
using std::string;
using std::vector;

using Acc = FixedString<30>;

static Block* db_block, *query_block;
static std::unordered_map<Acc, unsigned, Acc::Hash> acc2oid_db, acc2oid_query;

static Sequence get_query_seq(const string& acc) {
	return query_block->seqs()[acc2oid_query.at(acc)];
}

static Sequence get_db_seq(const string& acc) {
	return db_block->seqs()[acc2oid_db.at(acc)];
}

static SequenceSet get_seqs(const vector<string>& accs) {
	SequenceSet out;
	for (const string& s : accs) {
		try {
			out.reserve(get_db_seq(s).length());
		}
		catch (std::out_of_range&) {
			throw std::runtime_error("Target accession not found: " + s);
		}
	}
	out.finish_reserve();
	for (size_t i = 0; i < accs.size(); ++i) {
		Sequence s = get_db_seq(accs[i]);
		out.assign(i, s.data(), s.end());
	}
	return out;
}

static TextBuffer* view_query(const string& query_acc, const string& buf, SequenceFile& query_file, SequenceFile& target_file, Search::Config& cfg, Statistics& stats) {
	static const BlockId BATCH_SIZE = 1024;
	const vector<string> target_acc = Util::Tsv::extract_column(buf, 1);
	//SequenceSet targets = target_file.seqs_by_accession(target_acc.begin(), target_acc.end());
	//vector<Letter> query = query_file.seq_by_accession(query_acc);
	SequenceSet targets = get_seqs(target_acc);
	vector<Letter> query;
	try {
		query = get_query_seq(query_acc).copy();
	}
	catch (std::out_of_range&) {
		throw std::runtime_error("Query accession not found: " + query_acc);
	}
	if (cfg.query_masking != MaskingAlgo::NONE)
		Masking::get()(query.data(), query.size(), cfg.query_masking, 0);
	if (cfg.target_masking != MaskingAlgo::NONE)
		for (BlockId i = 0; i < targets.size(); ++i)
			Masking::get()(targets.ptr(i), targets.length(i), cfg.target_masking, 0);

	const auto query_comp = Stats::composition(Sequence(query));
	const int query_len = Stats::count_true_aa(Sequence(query));
	list<Hsp> hsp;

	for (BlockId t = 0; t < (BlockId)target_acc.size(); t += BATCH_SIZE) {
		const BlockId t1 = std::min(t + BATCH_SIZE, (BlockId)target_acc.size());

		vector<Stats::TargetMatrix> matrices;
		for (BlockId i = t; i < t1; ++i)
			matrices.emplace_back(query_comp, query_len, targets[i]);

		const HspValues v = HspValues::COORDS | HspValues::IDENT | HspValues::LENGTH;
		DP::Targets dp_targets;
		for (BlockId i = t; i < t1; ++i)
			if (targets.length(i) > 0)
				dp_targets[DP::BandedSwipe::bin(v, (Loc)query.size(), 0, 0, 0, 0, 0)].emplace_back(targets[i], targets[i].length(), i, &matrices[i - t]);

		DP::Params params{ Sequence(query), "", Frame(0), (Loc)query.size(), nullptr, DP::Flags::FULL_MATRIX, v, stats, nullptr };
		hsp.splice(hsp.end(), DP::BandedSwipe::swipe(dp_targets, params));
	}

	hsp.sort(Hsp::cmp_evalue);
	Blast_tab_format fmt;
	TranslatedSequence query_seq;
	TextBuffer* out = new TextBuffer();
	Output::Info info{ SeqInfo { 0, 0, query_acc.c_str(), "", query_seq.source().length(), query_seq.source(), Sequence()}, false, nullptr, *out, {} };
	for (Hsp& h : hsp) {
		h.query_source_range = h.query_range;
		fmt.print_match(HspContext(h,
			0,
			0,
			query_seq,
			query_acc.c_str(),
			0,
			0,
			target_acc[h.swipe_target].c_str(),
			0,
			0,
			Sequence()), info);
	}

	return out;
}

void view_tsv() {
	if (config.input_ref_file.size() > 1)
		throw std::runtime_error("Too many arguments for --in.");
	if (config.database.empty())
		throw std::runtime_error("Missing argument: database file (-d)");
	/*if (config.query_file.empty())
		throw std::runtime_error("Missing argument: query file (-q)");*/
	if (config.query_file.size() > 1)
		throw std::runtime_error("Too many arguments for query file (--query/-q)");

	task_timer timer("Opening the database file");
	unique_ptr<SequenceFile> db(SequenceFile::auto_create({ config.database }, SequenceFile::Flags::NO_FASTA));
	score_matrix = ScoreMatrix("blosum62", -1, -1, 1, 0);
	score_matrix.set_db_letters(config.db_size ? config.db_size : db->letters());
	Masking::instance = unique_ptr<Masking>(new Masking(score_matrix));

	timer.go("Opening the query file");
	unique_ptr<SequenceFile> query_file(SequenceFile::auto_create(config.query_file, SequenceFile::Flags::NO_FASTA));

	/*if (db->type() != SequenceFile::Type::BLAST) // || query_file->type() != SequenceFile::Type::BLAST)
		throw std::runtime_error("BLAST database required.");*/

	timer.go("Opening the input file");
	TextInputFile in(config.input_ref_file.front());

	timer.go("Opening the output file");
	OutputFile output_file(config.output_file);
	OutputWriter writer{ &output_file };
	output_sink.reset(new ReorderQueue<TextBuffer*, OutputWriter>(0, writer));

	timer.go("Loading database");
	db_block = db->load_seqs(SIZE_MAX, nullptr, SequenceFile::LoadFlags::ALL);

	timer.go("Loading queries");
	query_block = query_file->load_seqs(SIZE_MAX, nullptr, SequenceFile::LoadFlags::ALL);

	timer.go("Building accession mapping");
	unsigned n = db_block->ids().size();
	acc2oid_db.reserve(n);
	for (unsigned i = 0; i < n; ++i)
		acc2oid_db[Util::Seq::seqid(db_block->ids()[i], false)] = i;
	n = query_block->ids().size();
	acc2oid_query.reserve(n);
	for (unsigned i = 0; i < n; ++i)
		acc2oid_query[Util::Seq::seqid(query_block->ids()[i], false)] = i;

	timer.go("Computing alignments");
	size_t query_idx = 0;
	mutex mtx, mtx_out;
	Search::Config cfg;

	auto worker = [&] {
		try {
			string query;
			size_t q;
			Statistics stats;
			for (;;) {
				string buf;
				{
					lock_guard<mutex> lock(mtx);
					query = Util::Tsv::fetch_block(in, buf);
					q = query_idx++;
				}
				if (q % 1000 == 0)
					std::cout << "#Query = " << q << " time = " << timer.seconds() << endl;
				if (query.empty())
					return;
				TextBuffer* out = view_query(query, buf, *db, *db, cfg, stats);
				if (config.no_reorder) {
					{
						lock_guard<mutex> lock(mtx_out);
						output_file.write(out->data(), out->size());
					}
					delete out;
				}
				else
					output_sink->push(q, out);
			}
		}
		catch (std::exception& e) {
			exit_with_error(e);
		}
	};

	vector<thread> threads;
	for (int i = 0; i < config.threads_; ++i)
		threads.emplace_back(worker);
	for (auto& t : threads)
		t.join();

	timer.go("Closing the output file");
	output_file.close();
	output_sink.reset();
	delete db_block;
	delete query_block;
}