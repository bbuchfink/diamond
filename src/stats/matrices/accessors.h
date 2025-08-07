/****
DIAMOND protein aligner
Copyright (C) 2020 Max Planck Society for the Advancement of Science e.V.

Code developed by Benjamin Buchfink <benjamin.buchfink@tue.mpg.de>

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

#pragma once

namespace Stats {
    const StandardMatrix& getBlosum45();
    const StandardMatrix& getBlosum50();
    const StandardMatrix& getBlosum62();
    const StandardMatrix& getBlosum80();
    const StandardMatrix& getBlosum90();
    const StandardMatrix& getPam30();
    const StandardMatrix& getPam70();
    const StandardMatrix& getPam250();
}
