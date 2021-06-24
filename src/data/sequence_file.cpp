/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#include <iostream>
#include "sequence_file.h"
#include "../basic/masking.h"
#include "reference.h"
#include "dmnd/dmnd.h"
#include "../util/system/system.h"
#include "../util/util.h"
#include "../util/algo/partition.h"
#include "../util/sequence/sequence.h"
#include "../util/parallel/multiprocessing.h"
#ifdef WITH_BLASTDB
#include "blastdb/blastdb.h"
#endif

using std::cout;
using std::endl;
using std::setw;

const EMap<SequenceFile::Type> EnumTraits<SequenceFile::Type>::to_string = { {SequenceFile::Type::DMND, "Diamond database" }, {SequenceFile::Type::BLAST, "BLAST database"} };

static string dict_file_name(const size_t query_block, const size_t target_block) {
	const string file_name = append_label("ref_dict_", query_block) + append_label("_", target_block);
	return join_path(config.parallel_tmpdir, file_name);
}

void SequenceFile::load_block(size_t block_id_begin, size_t block_id_end, size_t pos, bool use_filter, const vector<uint64_t>* filtered_pos, bool load_ids, Block* block) {
	static const size_t MAX_LOAD_SIZE = 2 * GIGABYTES;
	seek_offset(pos);
	size_t load_size = 0;
	for (size_t i = block_id_begin; i < block_id_end; ++i) {
		bool seek = false;
		if (use_filter && (*filtered_pos)[i]) {
			pos = (*filtered_pos)[i];
			seek = true;
		}
		const size_t l = block->seqs_.length(i);
		load_size += l;
		read_seq_data(block->seqs_.ptr(i), l, pos, seek);
		if (load_ids)
			read_id_data(block->ids_.ptr(i), block->ids_.length(i));
		else
			skip_id_data();
		if (type_ == Type::DMND)
			Masking::get().remove_bit_mask(block->seqs_.ptr(i), block->seqs_.length(i));
		if (load_size > MAX_LOAD_SIZE) {
			close_weakly();
			reopen();
			load_size = 0;
		}
	}
}

size_t SequenceFile::dict_block(const size_t ref_block)
{
	return config.multiprocessing ? ref_block : 0;
}

Block* SequenceFile::load_seqs(const size_t max_letters, bool load_ids, const BitVector* filter, bool fetch_seqs, bool lazy_masking, const Chunk& chunk)
{
	task_timer timer("Loading reference sequences");
	reopen();

	if(max_letters == 0)
		seek_chunk(chunk);
	init_seqinfo_access();

	size_t database_id = tell_seq();
	size_t letters = 0, seqs = 0, id_letters = 0, seqs_processed = 0, filtered_seq_count = 0;
	vector<uint64_t> filtered_pos;
	Block* block = new Block(alphabet_);
	
	SeqInfo r = read_seqinfo();
	size_t offset = r.pos;
	bool last = false;
	if (type() == Type::BLAST && sequence_count() != sparse_sequence_count())
		filter = builtin_filter();
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
			if (fetch_seqs) {
				block->seqs_.reserve(r.seq_len);
			}
			if (load_ids) {
				const size_t id_len = this->id_len(r, r_next);
				id_letters += id_len;
				if (fetch_seqs)
					block->ids_.reserve(id_len);
			}
			++filtered_seq_count;
			block->block2oid_.push_back((unsigned)database_id);
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
		return block;

	if (fetch_seqs) {
		block->seqs_.finish_reserve();
		if (load_ids) block->ids_.finish_reserve();
		if (false && type_ == Type::BLAST && config.algo == Config::Algo::QUERY_INDEXED && config.threads_ > 1 && !use_filter) {
			assert(!use_filter);
			Partition<size_t> p(filtered_seq_count, config.threads_);
			vector<std::thread> t;
			for (size_t i = 0; i < p.parts; ++i)
				t.emplace_back(&SequenceFile::load_block, this, p.begin(i), p.end(i), offset + p.begin(i), use_filter, &filtered_pos, false, block);
			for (std::thread& i : t)
				i.join();
		}
		else
			load_block(0, filtered_seq_count, offset, use_filter, &filtered_pos, load_ids, block);
		
		timer.finish();
		block->seqs_.print_stats();
	}

	if (config.multiprocessing || config.global_ranking_targets)
		blocked_processing = true;
	else
		blocked_processing = seqs_processed < sequence_count();

	if (blocked_processing) // should be always
		close_weakly();

	if (lazy_masking)
		block->masked_.resize(filtered_seq_count, false);

	return block;
}

void SequenceFile::get_seq()
{
	std::map<string, string> seq_titles;
	if (!config.query_file.empty()) {
		TextInputFile list(config.single_query_file());
		while (list.getline(), !list.eof()) {
			const vector<string> t(tokenize(list.line.c_str(), "\t"));
			if (t.size() != 2)
				throw std::runtime_error("Query file format error.");
			seq_titles[t[0]] = t[1];
		}
		list.close();
	}

	vector<Letter> seq;
	string id;
	bool all = config.seq_no.size() == 0 && seq_titles.empty();
	std::set<size_t> seqs;
	if (!all)
		for (vector<string>::const_iterator i = config.seq_no.begin(); i != config.seq_no.end(); ++i)
			seqs.insert(atoi(i->c_str()) - 1);
	const size_t max_letters = config.chunk_size == 0.0 ? std::numeric_limits<size_t>::max() : (size_t)(config.chunk_size * 1e9);
	size_t letters = 0;
	TextBuffer buf;
	OutputFile out(config.output_file);
	for (size_t n = 0; n < sequence_count(); ++n) {
		read_seq(seq, id);
		std::map<string, string>::const_iterator mapped_title = seq_titles.find(Util::Seq::seqid(id.c_str(), false));
		if (all || seqs.find(n) != seqs.end() || mapped_title != seq_titles.end()) {
			buf << '>' << (mapped_title != seq_titles.end() ? mapped_title->second : id) << '\n';
			if (config.reverse) {
				Sequence(seq).print(buf, value_traits, Sequence::Reversed());
				buf << '\n';
			}
			else if (config.hardmasked) {
				Sequence(seq).print(buf, value_traits, Sequence::Hardmasked());
				buf << '\n';
			}
			else
				buf << Sequence(seq) << '\n';
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

SequenceFile::~SequenceFile()
{
	if (dict_file_) {
		dict_file_->close();
		dict_file_.reset();
	}
}

SequenceFile* SequenceFile::auto_create(Flags flags, Metadata metadata) {
	if (exists(config.database + ".pin") || exists(config.database + ".pal")) {
#ifdef WITH_BLASTDB
		if (config.multiprocessing)
			throw std::runtime_error("--multiprocessing is not compatible with BLAST databases.");
		if (config.target_indexed)
			throw std::runtime_error("--target-indexed is not compatible with BLAST databases.");
		return new BlastDB(config.database, metadata, flags);
#else
		throw std::runtime_error("This executable was not compiled with support for BLAST databases.");
#endif
	}
	config.database = auto_append_extension_if_exists(config.database, DatabaseFile::FILE_EXTENSION);
	if (DatabaseFile::is_diamond_db(config.database)) {
		return new DatabaseFile(config.database, metadata, flags);
	}
	else if (!flag_any(flags, Flags::NO_FASTA)) {
		message_stream << "Database file is not a DIAMOND or BLAST database, treating as FASTA." << std::endl;
		config.input_ref_file = { config.database };
		TempFile* db;
		DatabaseFile::make_db(&db);
		DatabaseFile* r(new DatabaseFile(*db));
		delete db;
		return r;
	}
	throw std::runtime_error("Database does not have a supported format.");
}

void SequenceFile::load_dict_block(InputFile* f, const size_t ref_block)
{
	while (load_dict_entry(*f, ref_block));
}

void SequenceFile::load_dictionary(const size_t query_block, const size_t ref_blocks)
{
	if (!dict_file_ && !config.multiprocessing)
		return;
	task_timer timer("Loading dictionary", 3);
	if (config.multiprocessing) {
		dict_oid_ = vector<vector<uint32_t>>(ref_blocks);
		reserve_dict(ref_blocks);
		for (size_t i = 0; i < ref_blocks; ++i) {
			InputFile f(dict_file_name(query_block, i));
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
		reserve_dict(0);
		InputFile f(*t);
		load_dict_block(&f, 0);
		if (dict_oid_.front().size() != next_dict_id_)
			throw std::runtime_error("Dictionary corrupted.");
		f.close_and_delete();
		dict_file_.reset();
	}
}

void SequenceFile::free_dictionary()
{
	dict_oid_.clear();
	dict_oid_.shrink_to_fit();
}

size_t SequenceFile::total_blocks() const {
	const size_t c = config.chunk_size * 1e9;
	return (this->letters() + c - 1) / c;
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
		block_to_dict_id_[block] = vector<uint32_t>(seq_count, DICT_EMPTY);
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

uint32_t SequenceFile::dict_id(size_t block, size_t block_id, size_t oid, size_t len, const char* id, const Letter* seq)
{
	auto it = block_to_dict_id_.find(block);
	if (it == block_to_dict_id_.end() || block_id >= it->second.size())
		throw std::runtime_error("Dictionary not initialized.");
	vector<uint32_t>& v = it->second;
	uint32_t n = v[block_id];
	if (n != DICT_EMPTY)
		return n;
	{
		std::lock_guard<std::mutex> lock(dict_mtx_);
		n = v[block_id];
		if (n != DICT_EMPTY)
			return n;
		n = next_dict_id_++;
		v[block_id] = n;
		write_dict_entry(block, oid, len, id, seq);
		return n;
	}
}

size_t SequenceFile::oid(uint32_t dict_id, const size_t ref_block) const
{
	const size_t b = dict_block(ref_block);
	if (b >= dict_oid_.size() || dict_id >= dict_oid_[b].size())
		throw std::runtime_error("Dictionary not loaded.");
	return dict_oid_[b][dict_id];
}

SequenceFile::SequenceFile(Type type, Alphabet alphabet, Flags flags):
	flags_(flags),
	type_(type),
	alphabet_(alphabet)
{}

void db_info() {
	if (config.database.empty())
		throw std::runtime_error("Missing option for database file: --db/-d.");
	SequenceFile* db = SequenceFile::auto_create(SequenceFile::Flags::NO_FASTA | SequenceFile::Flags::NO_COMPATIBILITY_CHECK);
	const std::streamsize w = 25;
	cout << setw(w) << "Database type  " << to_string(db->type()) << endl;
	cout << setw(w) << "Database format version  " << db->db_version() << endl;
	if(db->type() == SequenceFile::Type::DMND)
		cout << setw(w) << "Diamond build  " << db->program_build_version() << endl;
	cout << setw(w) << "Sequences  " << db->sequence_count() << endl;
	if (db->type() == SequenceFile::Type::BLAST && db->sequence_count() != db->sparse_sequence_count())
		cout << setw(w) << "Sequences (filtered) " << db->sparse_sequence_count() << endl;
	cout << setw(w) << "Letters  " << db->letters() << endl;
	db->close();
	delete db;
}