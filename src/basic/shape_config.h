/****
DIAMOND protein aligner
Copyright (C) 2013-2021 Max Planck Society for the Advancement of Science e.V.
                        Benjamin Buchfink
                        Eberhard Karls Universitaet Tuebingen
						
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
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <string>
#include "shape.h"

class ShapeConfig
{

public:
	
	ShapeConfig():
		n_ (0)
	{ }

	ShapeConfig(const std::vector<std::string>& codes, unsigned count):
		n_ (0)
	{
		unsigned maxShapes = count == 0 ? (unsigned)codes.size() : count;
		for (unsigned i = 0; i < maxShapes; ++i) {
			shapes_[n_] = Shape(codes[i].c_str(), i);
			if (shapes_[n_].weight_ != shapes_[0].weight_)
				throw std::runtime_error("Seed shape weight has to be uniform.");
			n_++;
		}
	}

	unsigned count() const
	{ return n_; }

	const Shape& operator[](size_t i) const
	{ return shapes_[i]; }

	friend std::ostream& operator<<(std::ostream&s, const ShapeConfig& cfg)
	{
		for (unsigned i = 0; i < cfg.n_; ++i)
			s << cfg.shapes_[i] << (i<cfg.n_-1?",":"");
		return s;
	}

	std::vector<uint32_t> patterns(unsigned begin, unsigned end) const {
		std::vector<uint32_t> v;
		for (unsigned i = begin; i < end; ++i)
			v.push_back(shapes_[i].mask_);
		return v;
	}

private:

	Shape shapes_[Const::max_shapes];
	unsigned n_;

};

extern ShapeConfig shapes;
extern unsigned shape_from, shape_to;