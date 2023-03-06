#include <array>
#include "file.h"

using std::array;
using std::runtime_error;

namespace Util { namespace Tsv {

template<typename It>
static Schema schema(It begin, const std::vector<FileColumn>& output_fields) {
	Schema schema;
	for (auto i = output_fields.cbegin(); i != output_fields.cend(); ++i)
		schema.push_back(begin[i->file][i->column]);
	return schema;
}

void join(File& file1, File& file2, int column1, int column2, const std::vector<FileColumn>& output_fields, File& out) {
	if (output_fields.empty())
		throw runtime_error("Join with empty output");
	array<File*, 2> files{ &file1, &file2 };
	array<int, 2> cols{ column1, column2 };
	array<Table, 2> tables{ file1.schema(), file2.schema() };
	for (int i = 0; i < 2; ++i)
		tables[i] = files[i]->read_record();
	while (!tables[0].empty() && !tables[1].empty()) {
		array<int64_t, 2> keys = { tables[0].front().get<int64_t>(column1), tables[1].front().get<int64_t>(column2) };
		if (keys[0] < keys[1])
			tables[0] = files[0]->read_record();
		else if (keys[1] < keys[0])
			tables[1] = files[1]->read_record();
		else {
			*out.out_file_ << tables[output_fields.front().file].front().get(output_fields.front().column);
			for (auto it = output_fields.begin() + 1; it < output_fields.end(); ++it)
				*out.out_file_ << '\t' << tables[it->file].front().get(it->column);
			*out.out_file_ << '\n';
			for (int i = 0; i < 2; ++i)
				tables[i] = files[i]->read_record();
		}
	}
}

File* join(File& file1, File& file2, int column1, int column2, const std::vector<FileColumn>& output_fields) {
	array<Schema, 2> schemas{ file1.schema_, file2.schema_ };
	File* out = new File(schema(schemas.begin(), output_fields), "",  Flags::TEMP);
	join(file1, file1, column1, column2, output_fields, *out);
	return out;
}

}}