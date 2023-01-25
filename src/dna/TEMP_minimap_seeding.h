#include <string>
#include <vector>

#include "TEMP_minimap_structures.h"
#include "../run/config.h"
#include "../data/block/block.h"
#include "../data/sequence_file.h"




std::vector<MinimizerHit> mainMap(const Search::Config &cfg,BlockId query_id,const Sequence &target);


