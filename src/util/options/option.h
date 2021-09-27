#pragma once

template<typename T>
struct Option : public T {
	Option():
		present_(false)
	{}
	bool present() const {
		return present_;
	}
	Option& operator=(const T& value) {
		*static_cast<T*>(this) = value;
		present_ = true;
		return *this;
	}
private:
	bool present_;
};