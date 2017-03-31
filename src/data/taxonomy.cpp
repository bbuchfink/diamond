/****
Copyright (c) 2017, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/

#include "taxonomy.h"
#include "../util/compressed_stream.h"
#include "../basic/config.h"
#include "../util/merge_sort.h"

Taxonomy taxonomy;

void Taxonomy::load()
{
	char acc[max_accesion_len + 2];
	unsigned taxid;
	Compressed_istream f(config.prot_accession2taxid);
	f.getline();
	
	while (!f.eof()) {
		f.getline();
		if (!f.line.empty() && sscanf(f.line.c_str(), "%*s%15s%u%*u", acc, &taxid) != 2) {
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