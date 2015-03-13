/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef MAKE_DB_H_
#define MAKE_DB_H_

#include <iostream>
#include "../basic/options.h"
#include "../data/reference.h"
#include "../basic/exceptions.h"
#include "../basic/statistics.h"
#include "../data/load_seqs.h"
#include "../util/seq_file_format.h"

template<class _val>
void make_db(_val)
{
	using std::cout;
	using std::endl;

	verbose_stream << "Database file = " << program_options::input_ref_file << endl;

	boost::timer::cpu_timer total;
	task_timer timer ("Opening the database file", true);
	Input_stream db_file (program_options::input_ref_file, true);
	timer.finish();

	ref_header.block_size = program_options::chunk_size;
	size_t chunk = 0;
	Output_stream main(program_options::database_file_name());
	main.write(&ref_header, 1);

	for(;;++chunk) {
		timer.go("Loading sequences");
		Sequence_set<Nucleotide>* ss;
		size_t n_seq = load_seqs<_val,_val,Single_strand>(db_file, FASTA_format<_val> (), ref_seqs<_val>::data_, ref_ids::data_, ss, (size_t)(program_options::chunk_size * 1e9));
		if(n_seq == 0)
			break;
		ref_header.letters += ref_seqs<_val>::data_->letters();
		ref_header.sequences += n_seq;
		const bool long_addressing = ref_seqs<_val>::data_->raw_len() > (size_t)std::numeric_limits<uint32_t>::max();
		ref_header.long_addressing = ref_header.long_addressing == true ? true : long_addressing;
		timer.finish();
		ref_seqs<_val>::data_->print_stats();

		timer.go("Building histograms");
		seed_histogram *hst = new seed_histogram (*ref_seqs<_val>::data_, _val());

		timer.go("Saving to disk");
		ref_seqs<_val>::data_->save(main);
		ref_ids::get().save(main);
		hst->save(main);

		timer.go("Deallocating sequences");
		delete ref_seqs<_val>::data_;
		delete ref_ids::data_;
		delete ss;
		delete hst;
	}

	timer.finish();
	ref_header.n_blocks = chunk;
	main.seekp(0);
	main.write(&ref_header, 1);
	main.close();

	verbose_stream << "Total time = " << boost::timer::format(total.elapsed(), 1, "%ws\n");
}

#endif /* MAKE_DB_H_ */
