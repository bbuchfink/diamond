/****
DIAMOND protein aligner
Copyright (C) 2020 QIAGEN A/S (Aarhus, Denmark)
Code developed by Patrick Ettenhuber <patrick.ettenhuber@qiagen.com>

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

#include "output_format.h"
#include "../data/queries.h"
#include "recursive_parser.h"

using std::string;

void Clustering_format::print_match(const HspContext& r, Output::Info& info)
{
	info.out.write((uint32_t) r.query_oid);
	info.out.write((uint32_t) r.subject_oid);
	RecursiveParser rp(&r, format.c_str());
	const double res = rp.evaluate();
	info.out.write(res);
}

Clustering_format::Clustering_format(const string* const format):
	OutputFormat(bin1, HspValues::NONE),
	format(RecursiveParser::clean_expression(format))
{	
	RecursiveParser rp(nullptr, format->c_str());
	rp.evaluate();
	const auto vars = rp.variables();
	for (const Variable* v : vars) {
		this->hsp_values |= v->hsp_values;
		this->flags |= v->flags;
	}
}
