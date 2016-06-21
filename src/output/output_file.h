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

#ifndef OUTPUT_FILE_H_
#define OUTPUT_FILE_H_

#include <string>
#include "output.h"
#include "output_buffer.h"

using std::string;

struct Block_output : public Buffered_file
{

	struct Iterator
	{
		unsigned block_;
		bool same_subject_;
		Intermediate_record info_;
		bool operator<(const Iterator &rhs) const
		{ return info_.query_id > rhs.info_.query_id ||
				(info_.query_id == rhs.info_.query_id && (rhs.same_subject_ ||
						(!rhs.same_subject_ && info_.score < rhs.info_.score))); }
	};

	bool next(Iterator &it, unsigned subject, unsigned query)
	{
		if(this->eof())
			return false;
		it.info_.read(*this);
		it.block_ = block_;
		it.same_subject_ = it.info_.subject_id == subject && it.info_.query_id == query;
		return true;
	}

	Block_output(unsigned ref_block, const Temp_file &tmp_file):
		Buffered_file (tmp_file),
		block_ (ref_block)
	{ }

private:

	const unsigned block_;

};

#endif /* OUTPUT_FILE_H_ */
