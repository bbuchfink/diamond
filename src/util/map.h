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

#ifndef MAP_H_
#define MAP_H_

template<typename _it, typename _key>
struct Map
{

	struct Iterator
	{
		Iterator(const _it& begin, const _it& parent_end):
			begin_ (begin),
			parent_end_ (parent_end),
			end_ (get_end())
		{ }
		void operator++()
		{ begin_ = end_; end_ = get_end(); }
		bool valid() const
		{ return begin_ < parent_end_; }
		_it& begin()
		{ return begin_; }
		_it& end()
		{ return end_; }
	private:
		_it get_end() const
		{
			_it i = begin_;
			while(i < parent_end_ && _key()(*i) == _key()(*begin_)) ++i;
			return i;
		}
		_it begin_, parent_end_, end_;
	};

	Map(const _it &begin, const _it &end):
		begin_ (begin),
		end_ (end)
	{ }

	Iterator begin()
	{ return Iterator(begin_, end_); }

private:
	_it begin_, end_;

};

#endif /* MAP_H_ */
