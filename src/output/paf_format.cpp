/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#include "output_format.h"
#include "util/sequence/sequence.h"
#include "stats/score_matrix.h"

void PAFFormat::print_query_intro(Output::Info& info) const
{
	if (info.unaligned) {
		info.out.write_until(info.query.title, Util::Seq::id_delimiters);
		info.out << "\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n";
	}
}

void PAFFormat::print_match(const HspContext& r, Output::Info& info)
{
	info.out.write_until(r.query_title.c_str(), Util::Seq::id_delimiters);
	info.out << '\t' << r.query.source().length() << '\t'
		<< r.query_source_range().begin_ << '\t'
		<< r.query_source_range().end_ - 1 << '\t'
		<< (Frame(r.frame()).strand == FORWARD ? '+' : '-') << '\t';
	print_title(info.out, r.target_title.c_str(), false, false, "<>");

	info.out << '\t' << r.subject_len << '\t'
		<< r.subject_range().begin_ << '\t'
		<< r.subject_range().end_ - 1 << '\t'
		<< r.identities() << '\t'
		<< r.length() << '\t';
#ifdef WITH_DNA
    if(config.command == ::Config::blastn){
        info.out << r.mapping_quality() << '\t'
        << "cm:i:" << r.n_anchors() << '\t';
    } else
        info.out << "255" << '\t';
#else
    info.out << "255" << '\t';
#endif
    info.out	<< "AS:i:" << (uint32_t)score_matrix.bitscore(r.score()) << '\t'
		<< "ZR:i:" << r.score() << '\t'
		<< "ZE:f:";
	info.out.print_e(r.evalue());
	info.out << '\n';
}