/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include <memory>
#include <string>
#include <numeric>
#include <limits>
#include "../util/binary_file.h"
#include "sorted_list.h"
#include "../basic/statistics.h"
#include "../data/seed_histogram.h"
#include "../util/hash_table.h"
#include "../util/hash_function.h"
#include "../basic/packed_loc.h"
#include "sequence_set.h"
#include "../util/ptr_vector.h"

using std::auto_ptr;

struct invalid_database_version_exception : public std::exception
{
	virtual const char* what() const throw()
	{
		return "Database was built with a different version of diamond as is incompatible.";
	}
};

struct Reference_header
{
	Reference_header() :
		unique_id(0x24af8a415ee186dllu),
		build(Const::build_version),
		db_version(current_db_version),
		sequences(0),
		letters(0)
	{ }
	uint64_t unique_id;
	uint32_t build, db_version;
	uint64_t sequences, letters, pos_array_offset;
#ifdef EXTRA
	Sequence_type sequence_type;
#endif
	enum { current_db_version = 1 };
};

extern Reference_header ref_header;

struct Database_format_exception : public std::exception
{
	virtual const char* what() const throw()
	{ return "Database file is not a DIAMOND database."; }
};

struct Database_file : public Input_stream
{
	Database_file():
		Input_stream (config.database)
	{
		if(this->read(&ref_header, 1) != 1)
			throw Database_format_exception ();
		if(ref_header.unique_id != Reference_header ().unique_id)
			throw Database_format_exception ();
		if(ref_header.build < min_build_required || ref_header.db_version != Reference_header::current_db_version)
			throw invalid_database_version_exception();
		if (ref_header.sequences == 0)
			throw std::runtime_error("Incomplete database file. Database building did not complete successfully.");
#ifdef EXTRA
		if(sequence_type(_val()) != ref_header.sequence_type)
			throw std::runtime_error("Database has incorrect sequence type for current alignment mode.");
#endif
		pos_array_offset = ref_header.pos_array_offset;
	}
	void rewind()
	{
		pos_array_offset = ref_header.pos_array_offset;
	}
	bool load_seqs();
	void get_seq();
	enum { min_build_required = 74 };

	size_t pos_array_offset;
};

void make_db();

struct ref_seqs
{
	static const Sequence_set& get()
	{ return *data_; }
	static Sequence_set& get_nc()
	{ return *data_; }
	static Sequence_set *data_;
};

struct ref_ids
{
	static const String_set<0>& get()
	{ return *data_; }
	static String_set<0> *data_;
};

extern Partitioned_histogram ref_hst;
extern unsigned current_ref_block;
extern bool blocked_processing;

inline size_t max_id_len(const String_set<0> &ids)
{
	size_t max (0);
	for(size_t i=0;i<ids.get_length(); ++i)
		max = std::max(max, find_first_of(ids[i].c_str(), Const::id_delimiters));
	return max;
}

string* get_allseqids(const char *s);

struct Ref_map
{
	Ref_map():
		next_ (0)
	{ }
	void init(unsigned ref_count)
	{
		const unsigned block = current_ref_block;
		if(data_.size() < block+1) {
			data_.resize(block+1);
			data_[block].insert(data_[block].end(), ref_count, std::numeric_limits<uint32_t>::max());
		}
	}
	uint32_t get(unsigned block, unsigned i)
	{
		uint32_t n = data_[block][i];
		if(n != std::numeric_limits<uint32_t>::max())
			return n;
		{
			mtx_.lock();
			n = data_[block][i];
			if (n != std::numeric_limits<uint32_t>::max()) {
				mtx_.unlock();
				return n;
			}
			n = next_++;
			data_[block][i] = n;
			len_.push_back((uint32_t)ref_seqs::get().length(i));
			const char *title = ref_ids::get()[i].c_str();
			if (config.salltitles)
				name_.push_back(new string(title));
			else if (config.sallseqid)
				name_.push_back(get_allseqids(title));
			else
				name_.push_back(get_str(title, Const::id_delimiters));
			mtx_.unlock();
		}
		return n;
	}

	unsigned length(uint32_t i) const
	{
		return len_[i];
	}
	const char* name(uint32_t i) const
	{
		return name_[i].c_str();
	}
	void init_rev_map()
	{
		rev_map_.resize(next_);
		unsigned n = 0;
		for (unsigned i = 0; i < data_.size(); ++i) {
			for (unsigned j = 0; j < data_[i].size(); ++j)
				if (data_[i][j] != std::numeric_limits<uint32_t>::max())
					rev_map_[data_[i][j]] = n + j;
			n += (unsigned)data_[i].size();
		}
	}
	unsigned original_id(unsigned i) const
	{
		return rev_map_[i];
	}
	uint32_t check_id(uint32_t i) const
	{
		if (i >= next_)
			throw std::runtime_error("Dictionary reference id out of bounds.");
		return i;
	}
private:
	tthread::mutex mtx_;
	vector<vector<uint32_t> > data_;
	vector<uint32_t> len_;
	Ptr_vector<string> name_;
	vector<uint32_t> rev_map_;
	uint32_t next_;
	friend void finish_daa(Output_stream&);
};

extern Ref_map ref_map;

#endif /* REFERENCE_H_ */
