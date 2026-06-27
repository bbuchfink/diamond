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
#include <stdexcept>
#include "cluster/multinode/len_sort.h"

int run_queue_stress_test();
int run_hit_buffer_stress_test();
void filestack();

namespace Test {

static void require(bool value, const char* message) {
	if (!value)
		throw std::runtime_error(message);
}

static void len_sort_block_limits() {
	using Cluster::Multinode::can_add_to_len_sorted_block;
	require(can_add_to_len_sorted_block(90, 9, 10, 100, 10, 1000), "Expected the last allowed sequence to fit.");
	require(!can_add_to_len_sorted_block(100, 10, 1, 1000, 10, 1000), "Expected sequence count cap to stop the block.");
	require(!can_add_to_len_sorted_block(90, 9, 11, 100, 10, 1000), "Expected letter cap to stop a non-empty block.");
	require(can_add_to_len_sorted_block(0, 0, 200, 100, 10, 1000), "Expected a single oversized sequence to form a block.");
	require(can_add_to_len_sorted_block(10, 2, 1, 1000, 10, 270), "Expected raw packed-position cap boundary to fit.");
	require(!can_add_to_len_sorted_block(10, 2, 2, 1000, 10, 270), "Expected raw packed-position cap to stop the block.");
	require(!can_add_to_len_sorted_block(0, 0, 800, 1000, 10, 1000), "Expected an unrepresentable sequence to fail.");
}

int run() {
	len_sort_block_limits();
	std::cerr << "Unit tests passed." << std::endl;
	return 0;
	//filestack();	
	
	//run_hit_buffer_stress_test();
	//run_queue_stress_test();	

}

}