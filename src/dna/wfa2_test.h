#include "../run/config.h"
#include "TEMP_minimap_structures.h"
#include "../basic/match.h"
#include "../align/extend.h"


namespace WaveExtension{
std::pair<std::vector<Extension::Match>, Extension::Stats> extend(const Search::Config &cfg,BlockId query_id);
}


