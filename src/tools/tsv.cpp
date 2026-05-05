/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#include "util/tsv/tsv.h"
#include "basic/config.h"
#include "util/log_stream.h"

using namespace Util::Tsv;
using std::endl;
using std::string;

void word_count() {
	const int64_t READ_SIZE = config.tsv_read_size;
	File f(Schema{ Type::STRING, Type::STRING }, config.query_file.front());
	verbose_stream << "File size: " << f.size() << " bytes" << endl;
	int64_t n = 0;
	for(;;) {
		const Table t = f.read(READ_SIZE, config.threads_);
		n += t.size();
		if (t.empty())
			break;
	}
	message_stream << n << endl;
}

void cut() {
	const int64_t READ_SIZE = config.tsv_read_size;
	Schema schema{ Type::STRING, Type::STRING };
	File f(schema, config.query_file.front());
	verbose_stream << "File size: " << f.size() << " bytes" << endl;
	auto callback = MapFunc([](const Record& r) -> Table {
		Table t(Schema{ Type::STRING });
		t.write_record(r.get<string>(0));
		return t;
		});
	File out(Schema{ Type::STRING }, "", Flags::WRITE);
	TaskTimer timer("", 3);
	for (;;) {
		timer.go("Loading data");
		const Table t = f.read(READ_SIZE, config.threads_);
		if (t.empty())
			break;
		//for (int64_t i = 0; i < t.size(); ++i)
			//std::cout << t[i].get(1) << std::endl;
		timer.go("Writing data");
		t.map(config.threads_, callback, out);
	}
}