#include <string>
#include <vector>
#pragma once

struct OptionsGroup;

struct OptionBase
{
	OptionBase(const std::string& id, char short_id, const std::string& desc, bool disabled, const OptionsGroup* group) :
		id(id),
		desc(desc),
		short_id(short_id),
		disabled(disabled),
		group(group)
	{}
	virtual void read(const std::vector<std::string>& v) = 0;
	virtual bool present() = 0;
	virtual void set_default() = 0;
	virtual ~OptionBase()
	{}
	const std::string id, desc;
	const char short_id;
	bool disabled;
	const OptionsGroup* group;
};

template<typename T>
struct Option : public T {
	Option():
		present_(false)
	{}
	bool present() const {
		return present_;
	}
	bool blank() const {
		return !present_;
	}
	Option& operator=(const T& value) {
		*static_cast<T*>(this) = value;
		present_ = true;
		return *this;
	}
	void require() {
		if(!present_)
			throw std::runtime_error("Missing parameter: --" + base_->id + "/-" + base_->short_id);
	}
	T get(const T& default_value) const {
		return present_ ? *this : default_value;
	}
	void unset() {
		present_ = false;
	}
private:
	bool present_;
	const OptionBase* base_;
	template<typename T2> friend void set_base_ptr(Option<T2>&, const OptionBase*);
};

template<>
struct Option<double> {
	Option() :
		present_(false)
	{}
	bool present() const {
		return present_;
	}
	bool blank() const {
		return !present_;
	}
	Option& operator=(const double value) {
		value_ = value;
		present_ = true;
		return *this;
	}
	void set_if_blank(const double value) {
		if (!present_)
			this->operator=(value);
	}
	void unset() {
		present_ = false;
	}
	operator double() const {
		if (!present_)
			throw std::runtime_error("Option::present");
		return value_;
	}
	double get(const double default_value) const {
		return present_ ? value_ : default_value;
	}
private:
	double value_;
	bool present_;
	const OptionBase* base_;
	template<typename T2> friend void set_base_ptr(Option<T2>&, const OptionBase*);
};

template<>
struct Option<int64_t> {
	Option() :
		present_(false)
	{}
	bool present() const {
		return present_;
	}
	bool blank() const {
		return !present_;
	}
	Option& operator=(const int64_t value) {
		value_ = value;
		present_ = true;
		return *this;
	}
	void set_if_blank(const int64_t value) {
		if (!present_)
			this->operator=(value);
	}
	void unset() {
		present_ = false;
	}
	operator int64_t() const {
		if (!present_)
			throw std::runtime_error("Option::present");
		return value_;
	}
	int64_t get(const int64_t default_value) const {
		return present_ ? value_ : default_value;
	}
private:
	int64_t value_;
	bool present_;
	const OptionBase* base_;
	template<typename T2> friend void set_base_ptr(Option<T2>&, const OptionBase*);
};