/****
DIAMOND protein aligner
Copyright (C) 2013-2018 Benjamin Buchfink <buchfink@gmail.com>

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

#include <utility>
#include "reference.h"
#include "ref_dictionary.h"
#include "../util/util.h"

using std::pair;

ReferenceDictionary ReferenceDictionary::instance_;

string* get_allseqids(const char *s)
{
	string *r = new string;
	const vector<string> t(tokenize(s, "\1"));
	for (vector<string>::const_iterator i = t.begin(); i != t.end(); ++i) {
		if (i != t.begin())
			r->append("\1");
		r->append(i->substr(0, find_first_of(i->c_str(), Const::id_delimiters)));
	}
	return r;
}

void ReferenceDictionary::clear()
{
	data_.clear();
	len_.clear();
	database_id_.clear();
	name_.clear();
	next_ = 0;
}

void ReferenceDictionary::init(unsigned ref_count, const vector<unsigned> &block_to_database_id)
{
	const unsigned block = current_ref_block;
	if (data_.size() < block + 1) {
		data_.resize(block + 1);
		data_[block].insert(data_[block].end(), ref_count, std::numeric_limits<uint32_t>::max());
	}
	block_to_database_id_ = &block_to_database_id;
}

uint32_t ReferenceDictionary::get(unsigned block, unsigned block_id)
{
	uint32_t n = data_[block][block_id];
	if (n != std::numeric_limits<uint32_t>::max())
		return n;
	{
		mtx_.lock();
		n = data_[block][block_id];
		if (n != std::numeric_limits<uint32_t>::max()) {
			mtx_.unlock();
			return n;
		}
		n = next_++;
		data_[block][block_id] = n;
		len_.push_back((uint32_t)ref_seqs::get().length(block_id));
		database_id_.push_back((*block_to_database_id_)[block_id]);
		const char *title = ref_ids::get()[block_id].c_str();
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

void ReferenceDictionary::build_lazy_dict(DatabaseFile &db_file)
{
	vector<bool> filter(db_file.ref_header.sequences);
	vector<pair<unsigned, unsigned> > m;
	const size_t dict_size = database_id_.size();
	m.reserve(dict_size);
	unsigned n = 0;
	for (vector<uint32_t>::const_iterator i = database_id_.begin(); i < database_id_.end(); ++i) {
		filter[*i] = true;
		m.push_back(std::make_pair(*i, n++));
	}
	db_file.rewind();
	vector<unsigned> block_to_database_id;
	db_file.load_seqs(block_to_database_id, std::numeric_limits<size_t>::max(), false, false, &filter);
	std::sort(m.begin(), m.end());
	dict_to_lazy_dict_id_.clear();
	dict_to_lazy_dict_id_.resize(dict_size);
	n = 0;
	for (vector<pair<unsigned, unsigned> >::const_iterator i = m.begin(); i < m.end(); ++i)
		dict_to_lazy_dict_id_[i->second] = n++;
}

/*void ReferenceDictionary::init_rev_map()
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
*/