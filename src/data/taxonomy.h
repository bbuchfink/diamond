/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include <vector>
#include <algorithm>
#include <string>
#include "../basic/const.h"
#include "../util/util.h"

using std::pair;
using std::string;

struct Taxonomy
{
	enum { max_accesion_len = 14 };
	struct Accession
	{
		Accession(const char *s)
		{
			strncpy(this->s, s, max_accesion_len);
		}
		Accession(const string &s)
		{
			string t(get_title(s));
			size_t i;
			if (t.compare(0, 6, "UniRef") == 0)
				t.erase(0, 9);
			else if ((i = t.find_first_of('|', 0)) != string::npos) {
				if (t.compare(0, 3, "gi|") == 0) {
					t.erase(0, t.find_first_of('|', i + 1) + 1);
					i = t.find_first_of('|', 0);
				}
				t.erase(0, i + 1);
				i = t.find_first_of('|', 0);
				if (i != string::npos)
					t.erase(i);
			}
			if (t.length() > max_accesion_len)
				this->s[0] = 0;
			else
				strncpy(this->s, t.c_str(), max_accesion_len);
		}
		bool operator<(const Accession &y) const
		{
			return strncmp(s, y.s, max_accesion_len) < 0;
		}
		bool match(const Accession &y) const
		{
			const void *p2 = memchr(y.s, '.', max_accesion_len);
			size_t n = max_accesion_len;
			if (p2 == 0) {
				const void *p1 = memchr(s, '.', max_accesion_len);
				if (p1)
					n = (const char*)p1 - s;
			}
			return strncmp(s, y.s, n) == 0;
		}
		friend std::ostream& operator<<(std::ostream &str, const Accession &x)
		{
			for (int i = 0; i < max_accesion_len && x.s[i] != 0; ++i)
				str << x.s[i];
			return str;
		}
		char s[max_accesion_len];
	};

	void load();

	unsigned get(const Accession &accession) const
	{
		std::vector<std::pair<Accession, unsigned> >::const_iterator i = std::lower_bound(accession2taxid_.begin(), accession2taxid_.end(), std::make_pair(accession, 0u));
		if (i < accession2taxid_.end() && i->first.match(accession))
			return i->second;
		else
			return 0;
	}

private:
	
	std::vector<std::pair<Accession, unsigned> > accession2taxid_;
};

extern Taxonomy taxonomy;