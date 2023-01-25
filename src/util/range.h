#pragma once

template<typename It>
struct Range {

	Range() {}

	Range(const It& begin, const It& end):
		begin_(begin),
		end_(end)
	{}

	It begin() const {
		return begin_;
	}

	It end() const {
		return end_;
	}

	size_t size() const {
		return end_ - begin_;
	}

private:

	It begin_, end_;

};