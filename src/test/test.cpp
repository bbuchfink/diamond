/****
DIAMOND protein sequence aligner
Copyright (C) 2012-2026 Benjamin J. Buchfink

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
// SPDX-License-Identifier: GPL-3.0-or-later

#include <iostream>

int run_queue_stress_test();
int run_hit_buffer_stress_test();
void filestack();

namespace Test {

int run() {
	std::cerr << "The test workflow is deprecated. Unit testing is available via CTest. No action was taken." << std::endl;
	return 0;
	//filestack();	
	
	//run_hit_buffer_stress_test();
	//run_queue_stress_test();	

}

}