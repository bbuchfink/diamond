/****
Copyright © 2013-2025 Benjamin J. Buchfink <buchfink@gmail.com>

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****/
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <vector>
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
		unsigned max_shapes = count == 0 ? (unsigned)codes.size() : std::min(count, (unsigned)codes.size());
		for (unsigned i = 0; i < max_shapes; ++i) {
			shapes_[n_] = Shape(codes[i].c_str(), i);
			if (shapes_[n_].weight_ != shapes_[0].weight_)
				throw std::runtime_error("Seed shape weight has to be uniform.");
			n_++;
		}
	}

	int count() const
	{ return n_; }

	const Shape& operator[](size_t i) const
	{ return shapes_[i]; }

	friend std::ostream& operator<<(std::ostream&s, const ShapeConfig& cfg)
	{
		for (unsigned i = 0; i < cfg.n_; ++i)
			s << cfg.shapes_[i] << (i < cfg.n_ - 1 ? "," : "");
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