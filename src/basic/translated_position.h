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

#ifndef TRANSLATED_POSITION_H_
#define TRANSLATED_POSITION_H_

struct TranslatedPosition
{

	TranslatedPosition(interval source, int translated, int frame):
		source(frame < 3 ? source.begin_ : source.end_ - 1),
		translated(translated),
		frame(frame)
	{}

	TranslatedPosition& operator++()
	{
		if (align_mode.query_translated) {
			if (frame < 3)
				source += 3;
			else
				source -= 3;
		}
		else
			++source;
		++translated;
		return *this;
	}

	int inner_end() const
	{
		if (align_mode.query_translated) {
			if (frame < 3)
				return source - 1;
			else
				return source + 1;
		}
		else
			return source - 1;
	}

	int source, translated, frame;

};

#endif