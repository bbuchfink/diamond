#include <string>
#include <limits.h>
#include <math.h>
#include <algorithm>
#include "edge_vec.h"
#include "../../util/io/text_input_file.h"
#include "../../util/string/tokenizer.h"
#include "../algo/merge_sort.h"
#include "../../basic/config.h"
#include "../io/temp_file.h"

using std::string;
using std::array;
using std::vector;

namespace Util { namespace Algo { namespace UPGMA_MC {

DistType dist_type() {
	return config.upgma_dist == "bitscore" ? DistType::BITSCORE : DistType::EVALUE;
}

EdgeVec::EdgeVec(const char *file):
	current_bucket_(-1),
	i_(0),
	size_(0),
	node_count_(0)
{
	array<vector<CompactEdge>, BUCKET_COUNT> buffers;

	const DistType dt = dist_type();

	if (config.upgma_input == "bin") {
		InputFile in(file, InputFile::NO_AUTODETECT);
		uint32_t query, target, max = 0;
		double score;
		try {
			for (;;) {
				in >> target;
				//std::cerr << "target=" << target << ' ' << std::numeric_limits<uint32_t>::max() << std::endl;
				if (target == std::numeric_limits<uint32_t>::max()) {
					in >> query;
					max = std::max(max, query);
					//std::cerr << "query=" << query << std::endl;
					continue;
				}
				max = std::max(max, target);
				in >> score;
				const double dist = -score;
				const int b = bucket(dist, dt);
				buffers[b].push_back({ (int)query, (int)target, dist });
				if (buffers[b].size() == 4096) {
					temp_files[b].write(buffers[b].data(), buffers[b].size());
					buffers[b].clear();
				}
				++size_;
			}
		}
		catch (EndOfStream &) {}
		node_count_ = max + 1;
		in.close();
	}
	else {
		TextInputFile in(file);
		string query, target;
		double evalue, bitscore;
		int qlen, slen;
		while (in.getline(), !in.eof()) {
			String::Tokenizer t(in.line, "\t");
			t >> query >> target;
			if (dt == DistType::BITSCORE)
				t >> bitscore >> qlen >> slen;
			else
				t >> evalue;
			auto i = acc2idx.emplace(query, (int)acc2idx.size()), j = acc2idx.emplace(target, (int)acc2idx.size());
			if (i.second)
				idx2acc[i.first->second] = i.first->first;
			if (j.second)
				idx2acc[j.first->second] = j.first->first;
			if (i.first->second < j.first->second) {
				const double dist = (dt == DistType::BITSCORE) ? -bitscore / std::max(qlen, slen) : evalue;
				const int b = bucket(dist, dt);
				//std::cout << b << std::endl;
				buffers[b].push_back({ i.first->second, j.first->second, dist });
				if (buffers[b].size() == 4096) {
					temp_files[b].write(buffers[b].data(), buffers[b].size());
					buffers[b].clear();
				}
				++size_;
			}
		}
		in.close();
	}
	for (int b = 0; b < BUCKET_COUNT; ++b)
		temp_files[b].write(buffers[b].data(), buffers[b].size());
}

CompactEdge EdgeVec::get() {
	while (i_ >= buffer_.size()) {
		++current_bucket_;
		if (current_bucket_ >= BUCKET_COUNT)
			return { 0, 0, 0.0 };
		const size_t n = temp_files[current_bucket_].tell() / sizeof(CompactEdge);
		InputFile in(temp_files[current_bucket_]);
		buffer_.clear();
		buffer_.resize(n);
		in.read(buffer_.data(), n);
		in.close_and_delete();
		i_ = 0;
		merge_sort(buffer_.begin(), buffer_.end(), config.threads_);
	}
	return buffer_[i_++];
}

std::string EdgeVec::print(int idx) const {
	return (size_t)idx < idx2acc.size() ? idx2acc.find(idx)->second : std::to_string(idx);
}

}}}