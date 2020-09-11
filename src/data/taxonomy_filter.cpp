#include "taxonomy.h"

using std::string;

TaxonomyFilter::TaxonomyFilter(const string &include, const string &exclude, const TaxonList &list, TaxonomyNodes &nodes):
	BitVector(list.size())
{
	if (!include.empty() && !exclude.empty())
		throw std::runtime_error("Options --taxonlist and --taxon-exclude are mutually exclusive.");
	const bool e = !exclude.empty();
	const std::set<unsigned> taxon_filter_list(parse_csv(e ? exclude : include));
	if (taxon_filter_list.empty())
		throw std::runtime_error("Option --taxonlist/--taxon-exclude used with empty list.");
	if (taxon_filter_list.find(1) != taxon_filter_list.end() || taxon_filter_list.find(0) != taxon_filter_list.end())
		throw std::runtime_error("Option --taxonlist/--taxon-exclude used with invalid argument (0 or 1).");
	for (size_t i = 0; i < list.size(); ++i)
		if (nodes.contained(list[i], taxon_filter_list) ^ e)
			set(i);
}