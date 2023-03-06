#include "../run/config.h"
#include "../basic/match.h"
#include "../align/extend.h"


namespace WaveExtension{
std::pair<std::vector<Extension::Match>, Extension::Stats> extend(const Search::Config &cfg,BlockId query_id);
}


