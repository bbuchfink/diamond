#ifndef RANGE_H_
#define RANGE_H_

template<typename _it>
struct Range {

	Range() {}

	Range(const _it &begin, const _it &end):
		begin_(begin),
		end_(end)
	{}

	_it begin() const {
		return begin_;
	}

	_it end() const {
		return end_;
	}

	size_t size() const {
		return end_ - begin_;
	}

private:

	_it begin_, end_;

};

#endif