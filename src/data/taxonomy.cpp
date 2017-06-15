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

#include <stdio.h>
#include "taxonomy.h"
#include "../util/compressed_stream.h"
#include "../basic/config.h"
#include "../util/merge_sort.h"

Taxonomy taxonomy;

string& get_accession(string &t)
{
	size_t i;
	if (t.compare(0, 6, "UniRef") == 0)
		t.erase(0, 9);
	else if ((i = t.find_first_of('|', 0)) != string::npos) {
		if (t.compare(0, 3, "gi|") == 0) {
			t.erase(0, t.find_first_of('|', i + 1) + 1);
			i = t.find_first_of('|', 0);
		}
		t.erase(0, i + 1);
		i = t.find_first_of('|', 0);
		if (i != string::npos)
			t.erase(i);
	}
	return t;
}

void Taxonomy::load()
{
	char acc[max_accesion_len + 2];
	unsigned taxid;
	Compressed_istream f(config.prot_accession2taxid);
	f.getline();
	
	while (!f.eof() && (f.getline(), !f.line.empty())) {
		if (sscanf(f.line.c_str(), "%*s%15s%u%*u", acc, &taxid) != 2) {
			//std::cout << f.line << endl;
			throw std::runtime_error("Invalid taxonomy mapping file format.");
		}
		if (strlen(acc) > max_accesion_len) {
			//std::cout << f.line << endl;
			throw std::runtime_error("Accession exceeds supported length.");
		}
		accession2taxid_.push_back(std::make_pair(Accession(acc), taxid));
		/*if (f.line_count % 10000 == 0)
			std::cout << f.line_count << endl;*/
	}
	merge_sort(accession2taxid_.begin(), accession2taxid_.end(), config.threads_);
}