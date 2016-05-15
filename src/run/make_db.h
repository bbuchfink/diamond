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

#include <limits>
#include <iostream>
#include "../basic/config.h"
#include "../data/reference.h"
#include "../basic/statistics.h"
#include "../data/load_seqs.h"
#include "../util/seq_file_format.h"

void make_db()
{
	using std::cout;
	using std::endl;

	message_stream << "Database file: " << config.input_ref_file << endl;
	message_stream << "Block size: " << (size_t)(config.chunk_size * 1e9) << endl;

	Timer total;
	total.start();
	task_timer timer ("Opening the database file", true);
	Compressed_istream db_file (config.input_ref_file);
	timer.finish();

	ref_header.block_size = config.chunk_size;
#ifdef EXTRA
	ref_header.sequence_type = sequence_type(_val ());
#endif
	size_t chunk = 0;
	Output_stream main(config.database, false);
	main.write(&ref_header, 1);

	for(;;++chunk) {
		timer.go("Loading sequences");
		Sequence_set* ss;
		size_t n_seq = load_seqs(db_file, FASTA_format(), (Sequence_set**)&ref_seqs::data_, ref_ids::data_, ss, (size_t)(config.chunk_size * 1e9));
		if(n_seq == 0)
			break;
		ref_header.letters += ref_seqs::data_->letters();
		ref_header.sequences += n_seq;
		const bool long_addressing = ref_seqs::data_->raw_len() > (size_t)std::numeric_limits<uint32_t>::max();
		ref_header.long_addressing = ref_header.long_addressing == true ? true : long_addressing;
		timer.finish();
		ref_seqs::data_->print_stats();

		timer.go("Building histograms");
		seed_histogram *hst = new seed_histogram (*ref_seqs::data_);

		timer.go("Saving to disk");
		ref_seqs::data_->save(main);
		ref_ids::get().save(main);
		hst->save(main);

		timer.go("Deallocating sequences");
		delete ref_seqs::data_;
		delete ref_ids::data_;
		delete ss;
		delete hst;
	}

	timer.finish();
	ref_header.n_blocks = (unsigned)chunk;
	main.seekp(0);
	main.write(&ref_header, 1);
	main.close();

	message_stream << "Total time = " << total.getElapsedTimeInSec() << "s" << std::endl;
}

#endif /* MAKE_DB_H_ */
