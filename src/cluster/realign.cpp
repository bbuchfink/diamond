/****
DIAMOND protein aligner
Copyright (C) 2022 Max Planck Society for the Advancement of Science e.V.

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

#include <float.h>
#include "../basic/config.h"
#include "../util/io/text_input_file.h"
#include "../util/log_stream.h"
#include "../data/sequence_file.h"
#include "cluster.h"
#include "../basic/match.h"
#include "../output/output_format.h"

using std::endl;
using std::vector;
using std::unique_ptr;
using std::function;
using std::runtime_error;
using std::string;

static const vector<string> DEFAULT_FORMAT = { "6", "qseqid", "sseqid", "approx_pident", "qstart", "qend", "sstart", "send", "evalue", "bitscore" };

namespace Cluster {

void realign() {
	config.database.require();
	config.clustering.require();

	OutputFile out(config.output_file);
	if (config.output_format.empty())
		config.output_format = DEFAULT_FORMAT;
	unique_ptr<OutputFormat> output_format(get_output_format());
	if (output_format->code != OutputFormat::blast_tab)
		throw runtime_error("The realign workflow only supports tabular output format.");
	if (Blast_tab_format::header_format(Config::cluster) == Header::SIMPLE)
		dynamic_cast<Blast_tab_format*>(output_format.get())->output_header(out, true);

	task_timer timer("Opening the database");
	unique_ptr<SequenceFile> db(SequenceFile::auto_create({ config.database }, SequenceFile::Flags::NEED_LETTER_COUNT | SequenceFile::Flags::ACC_TO_OID_MAPPING, SequenceFile::Metadata()));
	score_matrix.set_db_letters(config.db_size ? config.db_size : db->letters());
	config.max_evalue = DBL_MAX;
	timer.finish();
	message_stream << "#Database sequences: " << db->sequence_count() << ", #Letters: " << db->letters() << endl;
	if (flag_any(db->format_flags(), SequenceFile::FormatFlags::TITLES_LAZY))
		db->init_random_access(0, 0, false);

	FlatArray<OId> clusters;
	vector<OId> centroids;
	tie(clusters, centroids) = read<OId>(config.clustering, *db, CentroidSorted());
	message_stream << "Found " << centroids.size() << " centroids, " << clusters.data_size() << " mappings in input file." << endl;

	TextBuffer buf;
	function<void(const HspContext&)> format_output([&buf, &out, &output_format](const HspContext& h) {
		Output::Info info{ SeqInfo(), false, nullptr, buf, {} };
		info.query.title = h.query_title.c_str();
		output_format->print_match(h, info);
		out.write(buf.data(), buf.size());
		buf.clear();
	});
	realign(clusters, centroids, *db, format_output, output_format->hsp_values);

	timer.go("Freeing memory");
	db->close();
	out.close();
}

}