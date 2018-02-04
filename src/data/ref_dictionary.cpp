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

#include "reference.h"
#include "ref_dictionary.h"
#include "../util/util.h"

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

void ReferenceDictionary::init(unsigned ref_count)
{
	const unsigned block = current_ref_block;
	if (data_.size() < block + 1) {
		data_.resize(block + 1);
		data_[block].insert(data_[block].end(), ref_count, std::numeric_limits<uint32_t>::max());
	}
}

uint32_t ReferenceDictionary::get(unsigned block, unsigned i)
{
	uint32_t n = data_[block][i];
	if (n != std::numeric_limits<uint32_t>::max())
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

void ReferenceDictionary::init_rev_map()
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
