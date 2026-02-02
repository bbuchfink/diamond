/****
Copyright © 2012-2026 Benjamin J. Buchfink <buchfink@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
****/
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <iostream>
#include <thread>
#define _REENTRANT
#include "ips4o/ips4o.hpp"
#include "util/parallel/simple_thread_pool.h"
#include "blastdb/blastdb.h"
#include "sequence_file.h"
#include "masking/masking.h"
#include "dmnd/dmnd.h"
#include "util/system/system.h"
#include "util/sequence/sequence.h"
#include "util/parallel/multiprocessing.h"
#include "basic/config.h"
#include "fasta/fasta_file.h"
#include "util/string/tokenizer.h"
#include "util/tsv/tsv.h"
#include "util/log_stream.h"
#include "util/data_structures/queue.h"

using std::string;
using std::endl;
using std::setw;
using std::thread;
using std::pair;
using std::tie;
using std::runtime_error;
using std::ofstream;
using std::greater;
using std::numeric_limits;
using std::tuple;
using std::get;
using std::unique_ptr;
using std::function;
using std::vector;
using std::set;
using std::unordered_map;
using std::map;
using std::atomic;

static constexpr int64_t CHECK_FOR_DNA_COUNT = 10;
const char* const SequenceFile::SEQID_HDR = "seqid";
const DictId SequenceFile::DICT_EMPTY = std::numeric_limits<DictId>::max();

const EMap<SequenceFile::Type> EnumTraits<SequenceFile::Type>::to_string = {
	{SequenceFile::Type::DMND, "Diamond database" },
	{SequenceFile::Type::BLAST, "BLAST database"},
	{SequenceFile::Type::FASTA, "FASTA file"},
	{SequenceFile::Type::BLOCK, ""}
};

static string dict_file_name(const size_t query_block, const size_t target_block) {
	const string file_name = append_label("ref_dict_", query_block) + append_label("_", target_block);
	return join_path(config.parallel_tmpdir, file_name);
}

static size_t single_oid(const SequenceFile* f, const string& acc) {
	const vector<OId> oid = f->accession_to_oid(acc);
	if (oid.empty())
		throw AccessionNotFound();
	if (oid.size() > 1)
		throw std::runtime_error("Multiple oids for target accession: " + acc);
	return oid.front();
}

bool SequenceFile::files_synced() {
	throw OperationNotSupported();
}

void SequenceFile::init_write() {
	throw OperationNotSupported();
}

void SequenceFile::write_seq(const Sequence& seq, const std::string& id) {
	throw OperationNotSupported();
}

std::string SequenceFile::taxon_scientific_name(TaxId taxid) const {
	throw OperationNotSupported();
}

void SequenceFile::add_taxid_mapping(const vector<pair<OId, TaxId>>& taxids) {
	throw OperationNotSupported();
}

int SequenceFile::raw_chunk_no() const {
	throw OperationNotSupported();
}

pair<Block*, int64_t> SequenceFile::load_parallel(const uint64_t max_letters, const BitVector* filter, unordered_map<string, bool>* accs, const Chunk& chunk, bool load_taxids) {
	assert(chunk.n_seqs == 0);
	assert(config.threads_ > 0);
	const uint64_t letters = config.minichunk;
	const int t = std::min(config.load_threads, config.threads_), p = t >= 3 ? t - 2 : 1;
	const int raw_chunk_start = raw_chunk_no();
	int64_t seqs = 0;
	Block* block = new Block();
	Queue<RawChunk*> queue(p * 4, 1, p, nullptr);
	Queue<DecodedPackage*> output_queue(p * 4, p, 1, nullptr);
	const bool load_seqs = flag_any(flags(), SequenceFile::Flags::SEQS),
		load_titles = flag_any(flags(), SequenceFile::Flags::TITLES);
	SimpleThreadPool pool;
	auto worker = [&](const atomic<bool>& stop) {
		while (!stop) {
			RawChunk* c;
			if (!queue.wait_and_dequeue(c))
				break;
			DecodedPackage* pkg = c->decode(flags(), filter, accs);
			output_queue.enqueue(pkg);
			delete c;
		}
		output_queue.close();
		};
	auto writer = [&](const atomic<bool>& stop) {
		int next = raw_chunk_start;
		map<int, DecodedPackage*> backlog;
		while (!stop) {
			DecodedPackage* pkg;
			if (!output_queue.wait_and_dequeue(pkg))
				break;
			backlog.emplace(pkg->no, pkg);
			while (!backlog.empty() && backlog.begin()->first == next) {
				pkg = backlog.begin()->second;
				const size_t n_seqs = pkg->seqs.size();
				seqs += n_seqs;
				if (load_seqs) {
					block->seqs_.append(pkg->seqs);
				}
				if (load_titles) {
					block->ids_.append(pkg->ids);
				}
				if (load_taxids) {
					add_taxid_mapping(pkg->taxids);
				}
				block->block2oid_.insert(block->block2oid_.end(), pkg->oids.begin(), pkg->oids.end());
				delete pkg;
				backlog.erase(backlog.begin());
				++next;
			}
		}
		if (!backlog.empty())
			throw runtime_error("SequenceFile::load_parallel");
		};
	for (int i = 0; i < p; ++i)
		pool.spawn(worker);
	pool.spawn(writer);
	uint64_t block_letters = 0;
	RawChunk* rc;
	uint64_t bytes = 0;
	const Flags load_flags = flags() | (accs ? Flags::TITLES : Flags::NONE);
	do {
		rc = raw_chunk(std::min(letters, max_letters - block_letters), load_flags);
		if (rc->empty()) {
			delete rc;
			break;
		}
		block_letters += rc->letters();
		bytes += rc->bytes();
		queue.enqueue(rc);
	} while (!pool.stop() && block_letters < max_letters);
	queue.close();
	pool.join_all();
	if (seqs > 0) {
		if (load_seqs)
			block->seqs().finish_reserve();
		if (load_titles)
			block->ids().finish_reserve();
	}
	block->raw_bytes_ = bytes;
	return { block, seqs };
}

pair<Block*, int64_t> SequenceFile::load_twopass(const int64_t max_letters, const BitVector* filter, const Chunk& chunk) {
	init_seqinfo_access();

	OId database_id = tell_seq();
	int64_t letters = 0, seqs = 0, id_letters = 0, seqs_processed = 0, filtered_seq_count = 0;
	vector<int64_t> filtered_pos;
	Block* block = new Block();

	SeqInfo r = read_seqinfo();
	size_t offset = r.pos;
	bool last = false;
	const bool use_filter = filter && !filter->empty();

	auto goon = [&]() {
		if (max_letters > 0)
			return (r.seq_len > 0 && letters < max_letters);
		else
			return (seqs < chunk.n_seqs);
	};

	while (goon()) {
		SeqInfo r_next = read_seqinfo();
		if (!use_filter || filter->get(database_id)) {
			letters += r.seq_len;
			if (flag_any(flags(), SequenceFile::Flags::SEQS)) {
				block->seqs_.reserve(r.seq_len);
			}
			if (flag_any(flags(), SequenceFile::Flags::TITLES)) {
				const size_t id_len = this->id_len(r, r_next);
				id_letters += id_len;
				if (flag_any(flags(), SequenceFile::Flags::SEQS))
					block->ids_.reserve(id_len);
			}
			++filtered_seq_count;
			block->block2oid_.push_back(database_id);
			if (use_filter) {
				filtered_pos.push_back(last ? 0 : r.pos);
			}
			last = true;
		}
		else {
			last = false;
		}
		++database_id;
		++seqs_processed;
		r = r_next;
		++seqs;
	}

	putback_seqinfo();

	if (seqs == 0 || filtered_seq_count == 0)
		return { block, seqs_processed };
	const bool full_titles = bool(flags() & SequenceFile::Flags::FULL_TITLES);
	const bool all_seqids = bool(flags() & SequenceFile::Flags::ALL_SEQIDS);

	if (flag_any(flags(), SequenceFile::Flags::SEQS)) {
		block->seqs_.finish_reserve();
		if (flag_any(flags(), SequenceFile::Flags::TITLES)) block->ids_.finish_reserve();

		static const size_t MAX_LOAD_SIZE = 2 * GIGABYTES;
		if (use_filter && !flag_all(format_flags_, FormatFlags::SEEKABLE))
			throw OperationNotSupported();
		seek_offset(offset);
		for (BlockId i = 0; i < filtered_seq_count; ++i) {
			bool seek = false;
			if (use_filter && filtered_pos[i]) {
				offset = filtered_pos[i];
				seek = true;
			}
			const size_t l = block->seqs_.length(i);
			read_seq_data(block->seqs_.ptr(i), l, offset, seek);
			if (flag_any(flags(), SequenceFile::Flags::TITLES))
				read_id_data(block->block2oid_[i], block->ids_.ptr(i), block->ids_.length(i), all_seqids, full_titles);
			else
				skip_id_data();
			if (type_ == Type::DMND)
				Masking::get().remove_bit_mask(block->seqs_.ptr(i), block->seqs_.length(i));
		}
		
	}
	return { block, seqs_processed };
}

static int frame_mask() {
	if (config.query_strands == "both")
		return (1 << 6) - 1;
	else if (config.query_strands == "plus")
		return (1 << 3) - 1;
	else if (config.query_strands == "minus")
		return ((1 << 3) - 1) << 3;
	throw std::runtime_error("frame_mask");
}

std::pair<Block*, int64_t> SequenceFile::load_onepass(const int64_t max_letters, const BitVector* filter) {
	static const char* DNA_ERR = "The sequences are expected to be proteins but only contain DNA letters. Use the option --ignore-warnings to proceed.";
	vector<Letter> seq;
	string id;
	vector<char> qual;
	int64_t letters = 0, seq_count = 0;
	Block* block = new Block();
	const bool load_seqs = flag_any(flags(), SequenceFile::Flags::SEQS), load_titles = flag_any(flags(), SequenceFile::Flags::TITLES), preserve_dna = flag_any(flags(), SequenceFile::Flags::DNA_PRESERVATION);
	vector<char>* q = flag_any(flags(), SequenceFile::Flags::QUALITY) ? &qual : nullptr;
	OId oid = tell_seq();
	const bool first_block = oid == 0;
	const int frame_mask = ::frame_mask();
	const int64_t modulo = file_count();
	int looks_like_dna = 0;
	do {
		if (!read_seq(seq, id, q))
			break;
		if (seq.size() == 0)
			continue;
		if (filter && !filter->get(oid)) {
			++oid;
			continue;
		}

		letters += block->push_back(Sequence(seq), load_titles ? id.c_str() : nullptr, q, oid++, this->value_traits_.seq_type, frame_mask, !preserve_dna);

        ++seq_count;
		if (first_block && seq_count <= CHECK_FOR_DNA_COUNT && value_traits_.seq_type == SequenceType::amino_acid && Util::Seq::looks_like_dna(Sequence(seq)) && !config.ignore_warnings) {
			++looks_like_dna;
			if (looks_like_dna >= CHECK_FOR_DNA_COUNT)
				throw std::runtime_error(DNA_ERR);
		}
	} while (letters < max_letters || seq_count % modulo != 0);
	if (seq_count > 0 && looks_like_dna == seq_count)
		throw std::runtime_error(DNA_ERR);
	if (file_count() == 2 && !files_synced())
		throw std::runtime_error("Unequal number of sequences in paired read files.");
	block->seqs_.finish_reserve();
	if (value_traits_.seq_type == SequenceType::nucleotide)
		block->source_seqs_.finish_reserve();
	if (load_titles)
		block->ids_.finish_reserve();
	if (q)
		block->qual_.finish_reserve();
	return { block, seq_count };
}

size_t SequenceFile::dict_block(const size_t ref_block)
{
	return config.multiprocessing ? ref_block : 0;
}

Block* SequenceFile::load_seqs(const int64_t max_letters, const BitVector* filter, const Chunk& chunk)
{
	if (max_letters == 0)
		seek_chunk(chunk);

	Block* block;
	int64_t seqs_processed;
	if (type_ == Type::BLAST) {
		tie(block, seqs_processed) = load_parallel(max_letters, filter, nullptr, chunk, false);
	}
	else if (flag_any(format_flags_, FormatFlags::LENGTH_LOOKUP))
		tie(block, seqs_processed) = load_twopass(max_letters, filter, chunk);
	else {
		if (chunk.n_seqs)
			throw OperationNotSupported();
		tie(block, seqs_processed) = load_onepass(max_letters, filter);
	}

	if (block->empty())
		return block;

	if (flag_any(flags(), SequenceFile::Flags::LAZY_MASKING))
		block->masked_.resize(block->seqs_.size(), false);
	//block->seqs_.print_stats();
	return block;
}

void SequenceFile::get_seq()
{
	std::map<string, string> seq_titles;
	if (!config.query_file.empty()) {
		TextInputFile list(config.single_query_file());
		while (list.getline(), !list.eof()) {
			const vector<string> t(Util::String::tokenize(list.line.c_str(), "\t"));
			if (t.size() != 2)
				throw std::runtime_error("Query file format error.");
			seq_titles[t[0]] = t[1];
		}
		list.close();
	}

	vector<Letter> seq;
	string id;
	bool all = config.seq_no.size() == 0 && seq_titles.empty() && config.oid_list.empty();

	std::set<size_t> seqs;
	if (!all)
		for (vector<string>::const_iterator i = config.seq_no.begin(); i != config.seq_no.end(); ++i)
			seqs.insert(atoi(i->c_str()) - 1);
	if (!config.oid_list.empty()) {
		TextInputFile f(config.oid_list);
		OId oid;
		while (f.getline(), !f.line.empty() || !f.eof()) {
			Util::String::Tokenizer<Util::String::CharDelimiter>(f.line, Util::String::CharDelimiter('\t')) >> oid;
			seqs.insert(oid);
		}
		f.close();
	}
	if (!seqs.empty())
		message_stream << "#Selected sequences: " << seqs.size() << endl;

	const size_t max_letters = config.chunk_size == 0.0 ? std::numeric_limits<size_t>::max() : (size_t)(config.chunk_size * 1e9);
	size_t letters = 0;
	TextBuffer buf;
	OutputFile out(config.output_file);
	for (uint64_t n = 0; n < sequence_count(); ++n) {
		read_seq(seq, id);
		std::map<string, string>::const_iterator mapped_title = seq_titles.find(Util::Seq::seqid(id.c_str()));
		if (all || seqs.find(n) != seqs.end() || mapped_title != seq_titles.end()) {
			if (config.reverse) {
				buf << '>' << (mapped_title != seq_titles.end() ? mapped_title->second : id) << '\n';
				Sequence(seq).print(buf, value_traits, Sequence::Reversed());
				buf << '\n';
			}
			else if (config.hardmasked) {
				buf << '>' << (mapped_title != seq_titles.end() ? mapped_title->second : id) << '\n';
				Sequence(seq).print(buf, value_traits, Sequence::Hardmasked());
				buf << '\n';
			}
			else {
				Util::Seq::format(seq, id.c_str(), nullptr, buf, "fasta", amino_acid_traits, 80);
			}
		}
		out.write(buf.data(), buf.size());
		letters += seq.size();
		if (letters >= max_letters)
			break;
		seq.clear();
		id.clear();
		buf.clear();
	}

	out.close();
}

Util::Tsv::File* SequenceFile::make_seqid_list() {
	Util::Tsv::File* f = new Util::Tsv::File(Util::Tsv::Schema{ Util::Tsv::Type::STRING }, "", Util::Tsv::Flags::TEMP);
	vector<Letter> seq;
	string id;
	init_seq_access();
	for (uint64_t n = 0; n < sequence_count(); ++n) {
		read_seq(seq, id);
		f->write_record(Util::Seq::seqid(id.c_str()));
	}
	return f;
}

SequenceFile::~SequenceFile()
{
	if (dict_file_) {
		dict_file_->close();
		dict_file_.reset();
	}
}

static bool is_blast_db(const string& path) {
	if (exists(path + ".pin") || exists(path + ".pal") || ends_with(path, ".pal")) {
		if (config.multiprocessing)
			throw std::runtime_error("--multiprocessing is not compatible with BLAST databases.");
		if (config.target_indexed)
			throw std::runtime_error("--target-indexed is not compatible with BLAST databases.");
		return true;
	}
	return false;
}

SequenceFile* SequenceFile::auto_create(const vector<string>& path, Flags flags, const ValueTraits& value_traits) {
	if (path.size() == 1) {
		if (is_blast_db(path.front()))
			return new BlastDB(path.front(), flags, value_traits);
		const string a = auto_append_extension_if_exists(path.front(), DatabaseFile::FILE_EXTENSION);
		if (DatabaseFile::is_diamond_db(a))
			return new DatabaseFile(a, flags, value_traits);
	}
	if (!flag_any(flags, Flags::NO_FASTA)) {
		//message_stream << "Database file is not a DIAMOND or BLAST database, treating as FASTA." << std::endl;
		return new FastaFile(path, flags, value_traits);
	}
	throw std::runtime_error("Sequence file does not have a supported format.");
}

void SequenceFile::load_dict_block(InputFile* f, const size_t ref_block)
{
	while (load_dict_entry(*f, ref_block));
}

void SequenceFile::load_dictionary(const size_t query_block, const size_t ref_blocks)
{
	if (!dict_file_ && !config.multiprocessing)
		return;
	TaskTimer timer("Loading dictionary", 3);
	if (config.multiprocessing) {
		dict_oid_ = vector<vector<OId>>(ref_blocks);
		if (flag_any(flags_, Flags::SELF_ALN_SCORES))
			dict_self_aln_score_ = vector<vector<double>>(ref_blocks);
		reserve_dict(ref_blocks);
		for (size_t i = 0; i < ref_blocks; ++i) {
			InputFile f(dict_file_name(query_block, i), InputFile::NO_AUTODETECT);
			load_dict_block(&f, i);
			f.close_and_delete();
		}
	}
	else {
		TempFile* t = dynamic_cast<TempFile*>(dict_file_.get());
		if (!t)
			throw std::runtime_error("Failed to load dictionary file.");
		dict_oid_ = { {} };
		dict_oid_.front().reserve(next_dict_id_);
		if (flag_any(flags_, Flags::SELF_ALN_SCORES)) {
			dict_self_aln_score_ = { {} };
			dict_self_aln_score_.front().reserve(next_dict_id_);
		}
		reserve_dict(0);
		InputFile f(*t);
		load_dict_block(&f, 0);
		if ((DictId)dict_oid_.front().size() != next_dict_id_)
			throw std::runtime_error("Dictionary corrupted.");
		f.close_and_delete();
		dict_file_.reset();
	}
}

void SequenceFile::free_dictionary()
{
	dict_oid_.clear();
	dict_oid_.shrink_to_fit();
	dict_len_.clear();
	dict_len_.shrink_to_fit();
	dict_title_.clear();
	dict_title_.shrink_to_fit();
	dict_seq_.clear();
	dict_seq_.shrink_to_fit();
	dict_self_aln_score_.clear();
	dict_self_aln_score_.shrink_to_fit();
	block_to_dict_id_.clear();
}

size_t SequenceFile::total_blocks() const {
	const size_t c = (size_t)(config.chunk_size * 1e9);
	return (this->letters() + c - 1) / c;
}

SequenceSet SequenceFile::seqs_by_accession(const std::vector<std::string>::const_iterator begin, const std::vector<std::string>::const_iterator end)
{
	SequenceSet out;
	vector<size_t> oids;
	oids.reserve(end - begin);
	for (auto it = begin; it != end; ++it) {
		try {
			const size_t oid = single_oid(this, *it);
			oids.push_back(oid);
			out.reserve(seq_length(oid));
		}
		catch (AccessionNotFound&) {
			out.reserve(0);
			oids.push_back(SIZE_MAX);
		}		
	}
	out.finish_reserve();
	vector<Letter> seq;
	for (size_t i = 0; i < oids.size(); ++i) {
		if (oids[i] == SIZE_MAX)
			continue;
		seq_data(oids[i], seq);
		out.assign(i, seq.begin(), seq.end());
	}
	return out;
}

std::vector<Letter> SequenceFile::seq_by_accession(const std::string& acc)
{
	const size_t oid = single_oid(this, acc);
	vector<Letter> seq;
	seq_data(oid, seq);
	return seq;
}

void SequenceFile::init_dict(const size_t query_block, const size_t target_block)
{
	if (dict_file_)
		dict_file_->close();
	dict_file_.reset(config.multiprocessing ? new OutputFile(dict_file_name(query_block, target_block)) : new TempFile());
	next_dict_id_ = 0;
	dict_alloc_size_ = 0;
	block_to_dict_id_.clear();
}

void SequenceFile::init_dict_block(size_t block, size_t seq_count, bool persist)
{
	if(!persist)
		block_to_dict_id_.clear();
	if(block_to_dict_id_.find(block) == block_to_dict_id_.end())
		block_to_dict_id_[block] = vector<DictId>(seq_count, DICT_EMPTY);
}

void SequenceFile::close_dict_block(bool persist)
{
	if (config.multiprocessing) {
		dict_file_->close();
		dict_file_.reset();
	}
	if (!persist)
		block_to_dict_id_.clear();
}

DictId SequenceFile::dict_id(size_t block, size_t block_id, size_t oid, size_t len, const char* id, const Letter* seq, const double self_aln_score)
{
	auto it = block_to_dict_id_.find(block);
	if (it == block_to_dict_id_.end() || block_id >= it->second.size())
		throw std::runtime_error("Dictionary not initialized.");
	vector<DictId>& v = it->second;
	DictId n = v[block_id];
	if (n != DICT_EMPTY)
		return n;
	{
		std::lock_guard<std::mutex> lock(dict_mtx_);
		n = v[block_id];
		if (n != DICT_EMPTY)
			return n;
		n = next_dict_id_++;
		v[block_id] = n;
		write_dict_entry(block, oid, len, id, seq, self_aln_score);
		return n;
	}
}

size_t SequenceFile::oid(DictId dict_id, const size_t ref_block) const
{
	const size_t b = dict_block(ref_block);
	if (b >= dict_oid_.size() || dict_id >= (DictId)dict_oid_[b].size())
		throw std::runtime_error("Dictionary not loaded.");
	return dict_oid_[b][dict_id];
}

double SequenceFile::dict_self_aln_score(const size_t dict_id, const size_t ref_block) const
{
	const size_t b = dict_block(ref_block);
	if (b >= dict_self_aln_score_.size() || dict_id >= dict_self_aln_score_[b].size())
		throw std::runtime_error("Dictionary not loaded.");
	return dict_self_aln_score_[b][dict_id];
}

SequenceFile::SequenceFile(Type type, Flags flags, FormatFlags format_flags, const ValueTraits& value_traits):
	flags_(flags),
	format_flags_(format_flags),
	value_traits_(value_traits),
	type_(type)
{
	if (flag_any(flags_, Flags::OID_TO_ACC_MAPPING))
		//seqid_file_.reset(new File(Schema{ ::Type::INT64, ::Type::STRING }, "", ::Flags::TEMP)); // | ::Flags::RECORD_ID_COLUMN));
		seqid_file_.reset(new Util::Tsv::File(Util::Tsv::Schema{ Util::Tsv::Type::STRING }, "", Util::Tsv::Flags::TEMP));
}

void SequenceFile::write_dict_entry(size_t block, size_t oid, size_t len, const char* id, const Letter* seq, const double self_aln_score)
{
	OutputFile& f = *dict_file_;
	f.write((uint32_t)oid);
	if (flag_any(format_flags_, FormatFlags::DICT_LENGTHS))
		f.write((uint32_t)len);
	if (flag_any(format_flags_, FormatFlags::DICT_SEQIDS)) {
		f << id;
		dict_alloc_size_ += strlen(id);
	}
	if (flag_any(flags_, Flags::TARGET_SEQS))
		f.write(seq, len);
	if (flag_any(flags_, Flags::SELF_ALN_SCORES))
		f.write(self_aln_score);
}

bool SequenceFile::load_dict_entry(InputFile& f, size_t ref_block)
{
	const size_t b = dict_block(ref_block);
	uint32_t oid, len;
	string title;
	if(f.read(&oid, 1) == 0)
		return false;
	dict_oid_[b].push_back(oid);
	if (flag_any(format_flags_, FormatFlags::DICT_LENGTHS)) {
		f.read(len);
		dict_len_[b].push_back(len);
	}
	if (flag_any(format_flags_, FormatFlags::DICT_SEQIDS)) {
		f >> title;
		dict_title_[b].push_back(title.begin(), title.end());
	}
	if (flag_any(flags_, Flags::TARGET_SEQS)) {
		vector<Letter> v(len);
		f.read(v.data(), len);
		dict_seq_[b].push_back(v.begin(), v.end());
	}
	if (flag_any(flags_, Flags::SELF_ALN_SCORES)) {
		double self_aln_score;
		f.read(self_aln_score);
		dict_self_aln_score_[b].push_back(self_aln_score);
	}
	return true;
}

void SequenceFile::reserve_dict(const size_t ref_blocks)
{
	if (config.multiprocessing) {
		if (flag_any(format_flags_, FormatFlags::DICT_LENGTHS))
			dict_len_ = std::vector<std::vector<uint32_t>>(ref_blocks);
		if (flag_any(format_flags_, FormatFlags::DICT_SEQIDS))
			dict_title_ = std::vector<StringSet>(ref_blocks);
		if (flag_any(flags_, Flags::TARGET_SEQS))
			dict_seq_ = std::vector<SequenceSet>(ref_blocks);
	}
	else {
		if (flag_any(format_flags_, FormatFlags::DICT_LENGTHS)) {
			dict_len_ = { {} };
			dict_len_[0].reserve(next_dict_id_);
		}
		if (flag_any(format_flags_, FormatFlags::DICT_SEQIDS)) {
			dict_title_ = { {} };
			dict_title_[0].reserve(next_dict_id_, dict_alloc_size_);
		}
		if (flag_any(flags_, Flags::TARGET_SEQS)) {
			dict_seq_ = { {} };
			dict_seq_[0].reserve(next_dict_id_, 0);
		}
	}
}

std::string SequenceFile::dict_title(DictId dict_id, const size_t ref_block) const
{
	const size_t b = dict_block(ref_block);
	if (b >= dict_title_.size() || dict_id >= (DictId)dict_title_[b].size())
		throw std::runtime_error("Dictionary not loaded.");
	return dict_title_[b][dict_id];
}

Loc SequenceFile::dict_len(DictId dict_id, const size_t ref_block)
{
	const size_t b = dict_block(ref_block);
	if (b >= dict_len_.size() || dict_id >= (DictId)dict_len_[b].size())
		throw runtime_error("Dictionary not loaded.");
	return dict_len_[b][dict_id];
}

std::vector<Letter> SequenceFile::dict_seq(DictId dict_id, const size_t ref_block)
{
	const size_t b = dict_block(ref_block);
	if (b >= dict_seq_.size() || dict_id >= (DictId)dict_seq_[b].size())
		throw std::runtime_error("Dictionary not loaded.");
	Sequence s = dict_seq_[b][dict_id];
	return vector<Letter>(s.data(), s.end());
}

DbFilter* SequenceFile::filter_by_taxonomy(std::istream& filter, char delimiter, bool exclude)
{
	DbFilter *f = new DbFilter(sequence_count());
	set<TaxId> taxon_filter_list;
	string token;
	while (std::getline(filter, token, delimiter)) {
		taxon_filter_list.insert(std::stoi(token));
	}
	if (taxon_filter_list.empty())
		throw runtime_error("Option --taxonlist/--taxon-exclude used with empty list.");
	if (taxon_filter_list.find(1) != taxon_filter_list.end() || taxon_filter_list.find(0) != taxon_filter_list.end())
		throw runtime_error("Option --taxonlist/--taxon-exclude used with invalid argument (0 or 1).");
	for (OId i = 0; i < sequence_count(); ++i) {
		const bool c = contained(taxids(i), taxon_filter_list, exclude, exclude);
		if (c ^ exclude) {
			f->oid_filter.set(i);
			f->letter_count += seq_length(i);
		}
	}
	return f;
}

void SequenceFile::build_acc_to_oid() {
	acc2oid_.reserve(sequence_count());
	set_seqinfo_ptr(0);
	vector<Letter> seq;
	string id;
	for (OId i = 0; i < sequence_count(); ++i) {
		read_seq(seq, id);
		acc2oid_[Util::Seq::seqid(id.c_str())] = i;
	}
}

vector<OId> SequenceFile::accession_to_oid(const string& accession) const {
	try {
		return { acc2oid_.at(accession) };
	}
	catch (std::out_of_range&) {
		throw runtime_error("Accession not found in database: " + accession);
	}
}

void SequenceFile::init_random_access(const size_t query_block, const size_t ref_blocks, bool dictionary)
{
	if (dictionary)
		load_dictionary(query_block, ref_blocks);
}

std::string SequenceFile::seqid(OId oid, bool all, bool full_titles) {
	throw std::runtime_error("seqid");
	//if (oid >= acc_.size())
		//throw std::runtime_error("OId to accession mapping not available.");
	//return acc_[oid];
}

void SequenceFile::write_accession_list(const std::vector<bool>& oids, std::string& file_name) {
	ofstream f(file_name);
	const Util::Tsv::Table acc = seqid_file().read(config.threads_);
	for (OId i = 0; i < sequence_count(); ++i)
		if (!oids[i])
			f << acc[i].get<string>(0) << endl;
}

template<typename It>
std::vector<int64_t> SequenceFile::seq_offsets(It begin, It end) {
	assert(std::is_sorted(begin, end));
	vector<int64_t> r;
	if (end <= begin)
		return r;
	r.reserve(end - begin);
	set_seqinfo_ptr(0);
	init_seqinfo_access();
	const OId end_oid = *(end - 1) + 1;
	if (end_oid > sequence_count())
		throw runtime_error("OId out of bounds.");
	It it = begin;
	for (OId i = 0; i < end_oid; ++i) {
		SeqInfo info = read_seqinfo();
		if (i == *it) {
			if (it != begin && i - 1 == *(it - 1))
				r.push_back(-1);
			else
				r.push_back(info.pos);
			++it;
		}
	}
	return r;
}

template<typename It>
void SequenceFile::sub_db(It begin, It end, FastaFile* out) {
	vector<Letter> seq;
	string id;
	out->init_write();
	if (end <= begin)
		return;
	if (flag_any(format_flags_, FormatFlags::LENGTH_LOOKUP)) {
		const vector<int64_t> pos = seq_offsets(begin, end);
		for (int64_t p : pos) {
			if (p >= 0)
				seek_offset(p);
			read_seq(seq, id);
			out->write_seq(Sequence(seq), id);
		}
	}
	else {
		assert(std::is_sorted(begin, end));
		set_seqinfo_ptr(0);
		for (OId i = 0; i <= *(end - 1); ++i) {
			read_seq(seq, id);
			if (*begin == i) {
				out->write_seq(Sequence(seq), id);
				++begin;
			}
		}
	}
}

template void SequenceFile::sub_db<vector<SuperBlockId>::const_iterator>(vector<SuperBlockId>::const_iterator, vector<SuperBlockId>::const_iterator, FastaFile*);

template<typename It>
FastaFile* SequenceFile::sub_db(It begin, It end, const string& file_name) {
	FastaFile* f = new FastaFile(file_name, true, FastaFile::WriteAccess(), Flags::NEED_LENGTH_LOOKUP);
	sub_db(begin, end, f);
	return f;
}

template FastaFile* SequenceFile::sub_db<vector<OId>::const_iterator>(vector<OId>::const_iterator, vector<OId>::const_iterator, const string&);
template FastaFile* SequenceFile::sub_db<vector<SuperBlockId>::const_iterator>(vector<SuperBlockId>::const_iterator, vector<SuperBlockId>::const_iterator, const string&);

pair<int64_t, int64_t> SequenceFile::read_fai_file(const string& file_name, int64_t seqs, int64_t letters) {
	string acc;
	Loc len;
	TextInputFile fai(file_name);
	while (fai.getline(), !fai.line.empty() || !fai.eof()) {
		Util::String::Tokenizer<Util::String::CharDelimiter>(fai.line, Util::String::CharDelimiter('\t')) >> acc >> len;
		if (flag_any(flags_, Flags::ACC_TO_OID_MAPPING))
			acc2oid_[acc] = seqs;
		//if (flag_any(flags_, Flags::OID_TO_ACC_MAPPING))
//			acc_.push_back(acc.begin(), acc.end());
		++seqs;
		letters += len;
	}
	fai.close();
	return { seqs, letters };
}

void db_info() {
	using std::cout;
	if (config.database.empty())
		throw std::runtime_error("Missing option for database file: --db/-d.");
	SequenceFile* db = SequenceFile::auto_create({ config.database }, SequenceFile::Flags::NO_FASTA | SequenceFile::Flags::NO_COMPATIBILITY_CHECK);
	const std::streamsize w = 25;
	cout << setw(w) << "Database type  " << to_string(db->type()) << endl;
	cout << setw(w) << "Database format version  " << db->db_version() << endl;
	if(db->type() == SequenceFile::Type::DMND)
		cout << setw(w) << "Diamond build  " << db->program_build_version() << endl;
	cout << setw(w) << "Sequences  " << db->sequence_count() << endl;
	/*if (db->type() == SequenceFile::Type::BLAST && db->sequence_count() != db->sparse_sequence_count())
		cout << setw(w) << "Sequences (filtered) " << db->sparse_sequence_count() << endl;*/
	cout << setw(w) << "Letters  " << db->letters() << endl;
	db->close();
	delete db;
}

void SequenceFile::add_seqid_mapping(const std::string& id, OId oid) {
	const string acc = Util::Seq::seqid(id.c_str());
	if (flag_any(flags_, Flags::ACC_TO_OID_MAPPING)) {
		if (oid != (OId)acc2oid_.size())
			throw runtime_error("add_seqid_mapping");
		auto r = acc2oid_.emplace(acc, oid);
		if (!r.second)
			throw runtime_error("Accession is not unique in database file: " + acc);
	}
	if (flag_any(flags_, Flags::OID_TO_ACC_MAPPING)) {
		seqid_file_->write_record(acc);
	}
}

vector<tuple<FastaFile*, vector<OId>, Util::Tsv::File*>> SequenceFile::length_sort(int64_t block_size, function<int64_t(Loc)>& seq_size) {
	static const int64_t MIN_BLOCK_SIZE = 1;
	vector<tuple<FastaFile*, vector<OId>, Util::Tsv::File*>> files;
	init_seq_access();
	vector<Letter> seq;
	string id;
	vector<pair<Loc, OId>> lengths;
	lengths.reserve(sequence_count());
	for (OId i = 0; i < sequence_count(); ++i) {
		read_seq(seq, id);
		lengths.emplace_back((Loc)seq.size(), i);
	}
	ips4o::parallel::sort(lengths.begin(), lengths.end(), greater<pair<Loc, OId>>(), config.threads_);

	int64_t size = 0, seqs = 0, letters = 0;
	int block = 0;
	for (auto i = lengths.begin(); i != lengths.end(); ++i) {
		size += seq_size(i->first);
		letters += i->first;
		++seqs;
		i->first = block;
		if ((size >= block_size && letters >= MIN_BLOCK_SIZE) || seqs >= numeric_limits<SuperBlockId>::max()) {
			log_stream << "Super block " << block << " seqs=" << seqs << " letters=" << letters << endl;
			++block;
			size = 0;
			seqs = 0;
			letters = 0;
		}
	}
	if (size > 0) {
		log_stream << "Super block " << block << " seqs=" << seqs << " letters=" << letters << endl;
		++block;
	}
	ips4o::parallel::sort(lengths.begin(), lengths.end(), [](const pair<Loc, OId>& p1, const pair<Loc, OId>& p2) { return p1.second < p2.second; }, config.threads_);

	files.reserve(block);
	for (int i = 0; i < block; ++i) {
		//files.emplace_back(new FastaFile("", true, FastaFile::WriteAccess()), vector<OId>(),
		files.emplace_back(new FastaFile("", true, FastaFile::WriteAccess(), Flags::NEED_LENGTH_LOOKUP), vector<OId>(),
			new Util::Tsv::File(Util::Tsv::Schema{ Util::Tsv::Type::INT64 }, "", Util::Tsv::Flags::TEMP));
		get<0>(files.back())->init_write();
	}
	init_seq_access();
	//Util::Tsv::File map_file(Util::Tsv::Schema{ Util::Tsv::Type::INT64 }, "", Util::Tsv::Flags::WRITE_ACCESS);
	const string placeholer("X");
	for (OId i = 0; i < sequence_count(); ++i) {
		read_seq(seq, id);
		auto& f = files[lengths[i].first];
		get<0>(f)->write_seq(seq, placeholer);
		get<2>(f)->write_record(i);
		
	}
	for (auto f : files)
		get<0>(f)->set_seqinfo_ptr(0);
	return files;
}

Util::Tsv::File& SequenceFile::seqid_file() {
	seqid_file_->rewind();
	return *seqid_file_;
}

size_t SequenceFile::letters_filtered(const DbFilter& v) {
	size_t n = 0;
	for (OId i = 0; i < v.oid_filter.size(); ++i)
		if (v.oid_filter.get(i))
			n += seq_length(i);
	return n;
}

set<TaxId> SequenceFile::rank_taxid(const std::vector<TaxId>& taxid, int rank) {
	set<TaxId> r;
	for (TaxId i : taxid)
		r.insert(rank_taxid(i, rank));
	return r;
}

TaxId SequenceFile::get_parent(TaxId taxid) {
	throw OperationNotSupported();
}

int SequenceFile::rank(TaxId taxid) const {
	throw OperationNotSupported();
}

TaxId SequenceFile::rank_taxid(TaxId taxid, int rank) {
	static const int max = 64;
	int n = 0;
	while (true) {
		if (this->rank(taxid) == rank)
			return taxid;
		if (taxid <= 1)
			return 0;
		if (++n > max)
			throw runtime_error("Path in taxonomy too long (rank_taxid).");
		taxid = get_parent(taxid);
	}
	return 0;
}

vector<TaxId> SequenceFile::lineage(TaxId taxid) {
	vector<TaxId> out;
	int n = 0;
	while (true) {
		if (taxid <= 0)
			return {};
		if (taxid == 1)
			break;
		if (++n > MAX_LINEAGE)
			throw runtime_error("Path in taxonomy too long (TaxonomyNodes::lineage).");
		out.push_back(taxid);
		taxid = get_parent(taxid);
	}
	reverse(out.begin(), out.end());
	return out;
}

TaxId SequenceFile::get_lca(TaxId t1, TaxId t2)
{
	if (t1 == t2 || t2 <= 0)
		return t1;
	if (t1 <= 0)
		return t2;
	TaxId p = t2;
	set<TaxId> l;
	l.insert(p);
	int n = 0;
	do {
		p = get_parent(p);
		if (p <= 0)
			return t1;
		l.insert(p);
		if (++n > MAX_LINEAGE)
			throw runtime_error("Path in taxonomy too long (get_lca).");
	} while (p != t1 && p != 1);
	if (p == t1)
		return p;
	p = t1;
	n = 0;
	while (l.find(p) == l.end()) {
		p = get_parent(p);
		if (p <= 0)
			return t2;
		if (++n > MAX_LINEAGE)
			throw runtime_error("Path in taxonomy too long (get_lca).");
	}
	return p;
}

TaxId SequenceFile::max_taxid() const {
	throw OperationNotSupported();
}

bool SequenceFile::contained(TaxId query, const set<TaxId>& filter, bool include_invalid)
{
	static const int max = 64;
	if (get_parent(query) < 0)
		return include_invalid;
	if (cached_[query])
		return contained_[query];
	if (filter.find(1) != filter.end())
		return true;
	int n = 0;
	TaxId p = query;
	while (p > 1 && filter.find(p) == filter.end()) {
		p = get_parent(p);
		if (p <= 0)
			return include_invalid;
		if (++n > max)
			throw runtime_error("Path in taxonomy too long (contained).");
	}
	const bool contained = p > 1;
	TaxId q = query;
	while (set_cached(q, contained), q != p)
		q = get_parent(q);
	return contained;
}

bool SequenceFile::contained(const vector<TaxId>& query, const set<TaxId>& filter, bool all, bool include_invalid)
{
	static const int max = 64;
	if (filter.find(1) != filter.end())
		return true;
	for (vector<TaxId>::const_iterator i = query.begin(); i != query.end(); ++i) {
		const bool c = contained(*i, filter, include_invalid);
		if (c && !all)
			return true;
		if (!c && all)
			return false;
	}
	return all;
}

void SequenceFile::init_cache() {
	const TaxId m = max_taxid();
	cached_.insert(cached_.end(), m + 1, false);
	contained_.insert(contained_.end(), m + 1, false);
}

void SequenceFile::print_info() const {
	message_stream << "Database: " << config.database << ' ';
	message_stream << "(type: " << to_string(type()) << ", ";
	message_stream << "sequences: " << sequence_count() << ", ";
	message_stream << "letters: " << letters() << ')' << endl;
}

void SequenceFile::init_taxon_output_fields() {
	using namespace std::placeholders;
	auto callback = [](const TabularFormat&, const HspContext& r, Output::Info& info, Rank rank) {
		set<TaxId> tax_id = info.db->rank_taxid(info.db->taxids(r.subject_oid), rank);
		tax_id.erase(-1);
		tax_id.erase(0);
		tax_id.erase(1);
		if (tax_id.empty())
			info.out << "N/A";
		else
			print_taxon_names(tax_id.begin(), tax_id.end(), *info.db, info.out);
		};
	auto star = [](const TabularFormat&, Output::Info& info) {
		info.out << '*';
		};
	for (int i = 1; i < Rank::count; ++i) {
		const FieldId next = FieldId((int)TabularFormat::field_def.rbegin()->first + 1);
		OutputField f;
		f.description = "Unique subject ";
		f.description += Rank::names[i];
		f.description += "(s), separated by a ';'";
		f.flags = Output::Flags::IS_ARRAY | Output::Flags::NO_REALIGN;
		f.hsp_values = HspValues::NONE;
		f.id = next;
		f.key = "s";
		f.key += Rank::names[i];
		std::replace(f.key.begin(), f.key.end(), ' ', '_');
		if (!TabularFormat::field_def.emplace(next, f).second)
			throw runtime_error("Duplicate taxonomic rank.");
		TabularFormat::field_callbacks[next].match = std::bind(callback, _1, _2, _3, Rank(i));
		TabularFormat::field_callbacks[next].query_intro = star;
	}
}

RawChunk* SequenceFile::raw_chunk(size_t letters, SequenceFile::Flags flags) {
	throw OperationNotSupported();
}