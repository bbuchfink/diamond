/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
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
#include <memory>
#include <vector>
#include <map>
#include "../basic/statistics.h"
#include "legacy/query_mapper.h"
#include "../run/workflow.h"
#include "../search/hit.h"

struct Output_writer
{
	Output_writer(OutputFile* f) :
		f_(f)
	{ }
	void operator()(TextBuffer &buf)
	{
		f_->write(buf.data(), buf.size());
		buf.clear();
	}
private:
	OutputFile* const f_;
};

void align_queries(Consumer* output_file, Search::Config &cfg);
