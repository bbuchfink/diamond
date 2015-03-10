/****
Copyright (c) 2014, University of Tuebingen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****
Author: Benjamin Buchfink
****/

#ifndef SHAPE_CONFIG_H_
#define SHAPE_CONFIG_H_

#include "shape.h"

template<typename _val>
struct shape_codes
{
	static const char* str[Const::index_modes][Const::max_shapes];
};

template<> const char* shape_codes<Amino_acid>::str[][Const::max_shapes] = {
		{ "111101011101111", "111011001100101111", "1111001001010001001111", "111100101000010010010111", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },				// 4x12
		{ "1111011111",		// 16x9
		"111001101111",
		"11101100101011",
		"11010010111011",
		"111010100001111",
		"1110100011001011",
		"11100010100101011",
		"11011000001100111",
		"1101010010000010111",
		"11100001000100100111",
		"110110000100010001101",
		"1110000100001000101011",
		"1101010000010001001011",
		"1101001001000010000111",
		"1101000100100000100000111",
		"1110001000100000001010011" }
};

template<> const char* shape_codes<Nucleotide>::str[][Const::max_shapes] = {
		{ "111110111011110110111111", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
		{ "111011011011101011111", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
};

class shape_config
{

public:

	static shape_config instance;

	shape_config():
		n_ (0),
		mode_ (0)
	{ }

	template<class _val>
	shape_config(unsigned mode, const _val&):
		n_ (0),
		mode_ (mode-1)
	{
		unsigned maxShapes = program_options::shapes == 0 ? Const::max_shapes : program_options::shapes;
		for(unsigned i=0;i<maxShapes;++i)
			if(shape_codes<_val>::str[mode_][i])
				shapes_[n_++] = shape (shape_codes<_val>::str[mode_][i], i);
	}

	unsigned count() const
	{ return n_; }

	const shape& get_shape(unsigned i) const
	{ return shapes_[i]; }

	unsigned mode() const
	{ return mode_; }

	static const shape_config& get()
	{ return instance; }

private:

	shape shapes_[Const::max_shapes];
	unsigned n_, mode_;

};

shape_config shape_config::instance;

#endif /* SHAPE_CONFIG_H_ */
