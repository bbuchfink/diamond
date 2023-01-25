#include "../run/config.h"
#include "../data/seed_array.h"


namespace Dna{
class Index{
public:
    Index(Search::Config& cfg,char *ref_buffer);

private:
    std::unique_ptr<SeedArray> ref_idx_;
};
}