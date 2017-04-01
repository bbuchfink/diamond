/****
Copyright (c) 2017, Benjamin Buchfink
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
			if (t.compare(0, 6, "UniRef")) {
				t.erase(0, 9);
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
		char s[max_accesion_len];
	};

	void load();

	unsigned get(const Accession &accession) const
	{
		std::vector<std::pair<Accession, unsigned> >::const_iterator i = std::lower_bound(accession2taxid_.begin(), accession2taxid_.end(), std::make_pair(accession, 0u));
		if (i->first.match(accession))
			return i->second;
		else
			return 0;
	}

private:
	
	std::vector<std::pair<Accession, unsigned> > accession2taxid_;
};

extern Taxonomy taxonomy;