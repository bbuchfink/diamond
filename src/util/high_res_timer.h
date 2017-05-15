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

#include <stdint.h>
#ifdef _MSC_VER
#else
#include <time.h>
#endif

struct High_res_timer {
	
	High_res_timer():
#ifdef _MSC_VER
	time_(__rdtsc())
#else
	time_(0)
#endif
	{
	}

	uint64_t get() const
	{
#ifdef _MSC_VER
		return __rdtsc() - time_;
#else
		return 0;
#endif
	}

	uint64_t nanoseconds() {
#ifdef _MSC_VER
		return 0;		
#else
		return 0;
#endif
	}

	double microseconds()
	{
		return nanoseconds() / 1000.0;
	}

private:
	unsigned long long time_;

};