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
#include <vector>
#include <atomic>

template<typename _t>
struct DynamicIterator {

	DynamicIterator(size_t count):
		count(count)
	{}

	virtual _t operator++(int) = 0;
	virtual _t operator[](size_t i) = 0;
	virtual ~DynamicIterator() {}

	const size_t count;

};

template<typename _t>
struct VectorIterator : public DynamicIterator<_t> {

	VectorIterator(typename std::vector<_t>::const_iterator begin, typename std::vector<_t>::const_iterator end):
		DynamicIterator<_t>(end - begin),
		i_(0),
		begin_(begin)
	{}

	virtual _t operator++(int) override {
		const size_t j = i_++;
		if (j < DynamicIterator<_t>::count)
			return begin_[j];
		else
			return _t();
	}

	virtual _t operator[](size_t i) override {
		return begin_[i];
	}

	virtual ~VectorIterator() {}
	
private:

	std::atomic<size_t> i_;
	const typename std::vector<_t>::const_iterator begin_;

};

template<typename _t, typename _cont>
struct ContainerIterator : public DynamicIterator<_t> {

	ContainerIterator(_cont& container, size_t size) :
		DynamicIterator<_t>(size),
		container_(container),
		i_(0)
	{}

	virtual _t operator++(int) override {
		const size_t j = i_++;
		if (j < DynamicIterator<_t>::count)
			return _t(container_[j], j);
		else
			return _t();
	}

	virtual _t operator[](size_t i) override {
		return _t(container_[i], i);
	}

	virtual ~ContainerIterator() {}

private:

	_cont& container_;
	std::atomic<size_t> i_;

};