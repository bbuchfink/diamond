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
		{ return begin_ != parent_end_; }
		_it& begin()
		{ return begin_; }
		_it& end()
		{ return end_; }
	private:
		_it get_end() const
		{
			_it i = begin_;
			while(i != parent_end_ && _key()(*i) == _key()(*begin_)) ++i;
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
