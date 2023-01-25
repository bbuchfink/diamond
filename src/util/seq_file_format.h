/****
DIAMOND protein aligner
Copyright (C) 2013-2020 Max Planck Society for the Advancement of Science e.V.
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
#include <string>
#include <memory>
#include "../basic/value.h"
#include "io/text_input_file.h"

struct file_format_exception : public std::exception
{
	virtual const char* what() const throw()
	{
		return "Invalid input file format";
	}
};

struct Sequence_file_format
{

	virtual bool get_seq(std::string &id, std::vector<Letter> &seq, TextInputFile &s, const ValueTraits& value_traits, std::vector<char> *qual = nullptr) const = 0;
	virtual ~Sequence_file_format()
	{ }
	
};

struct FASTA_format : public Sequence_file_format
{

	FASTA_format()
	{ }

	virtual bool get_seq(std::string &id, std::vector<Letter> &seq, TextInputFile &s, const ValueTraits& value_traits, std::vector<char> *qual = nullptr) const override;

	virtual ~FASTA_format()
	{ }

};

struct FASTQ_format : public Sequence_file_format
{

	FASTQ_format()
	{ }

	virtual bool get_seq(std::string &id, std::vector<Letter> &seq, TextInputFile &s, const ValueTraits& value_traits, std::vector<char> *qual = nullptr) const override;

	virtual ~FASTQ_format()
	{ }

};

std::unique_ptr<const Sequence_file_format> guess_format(TextInputFile &file);