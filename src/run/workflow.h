/****
DIAMOND protein aligner
Copyright (C) 2016-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
						
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
#include "../basic/config.h"
#include "../util/data_structures/bit_vector.h"

struct DatabaseFile;
struct Consumer;
struct TextInputFile;

namespace Workflow { 
namespace Search {

struct Options {
	Options():
		self(config.self),
		db(nullptr),
		consumer(nullptr),
		query_file(nullptr),
		db_filter(nullptr)
	{}
	bool self;
	DatabaseFile *db;
	Consumer *consumer;
	TextInputFile *query_file;
	const BitVector* db_filter;
};

void run(const Options &options);

}
}
