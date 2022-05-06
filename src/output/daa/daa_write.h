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

#pragma once
#include "../output.h"
#include "daa_file.h"
#include "../data/sequence_file.h"

void init_daa(OutputFile& f);

size_t write_daa_query_record(TextBuffer& buf, const char* query_name, const Sequence& query);

void finish_daa_query_record(TextBuffer& buf, size_t seek_pos);

void write_daa_record(TextBuffer& buf, const IntermediateRecord& r);

void write_daa_record(TextBuffer& buf, const Hsp& match, uint32_t subject_id);

void finish_daa(OutputFile& f, const SequenceFile& db);

void finish_daa(OutputFile& f, DAA_file& daa_in);

void finish_daa(OutputFile& f, DAA_file& daa_in, const StringSet& seq_ids, const std::vector<uint32_t>& seq_lens, int64_t query_count);

