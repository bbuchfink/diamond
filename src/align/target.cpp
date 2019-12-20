/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen

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

#include <algorithm>
#include <vector>
#include "../basic/config.h"
#include "target.h"
#include "../basic/diagonal_segment.h"
#include "../dp/ungapped.h"

using std::vector;

namespace Extension {

	/*void Target::inner_culling()
	{
		hsps.sort();
		if (hsps.size() > 0)
			filter_score = hsps.front().score;
		else
			filter_score = 0;
		for (list<Hsp>::iterator i = hsps.begin(); i != hsps.end();) {
			if (i->is_enveloped_by(hsps.begin(), i, 0.5) || (int)i->score < cutoff)
				i = hsps.erase(i);
			else
				++i;
		}

	}*/

}