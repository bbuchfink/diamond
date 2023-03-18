#include "../util/tsv/file.h"
#include "../basic/config.h"
#include "../util/log_stream.h"

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