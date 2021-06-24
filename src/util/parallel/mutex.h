#pragma once
#include <mutex>

struct Sync {};
struct Async {};

template<typename Tag>
struct Mutex {};

template<>
struct Mutex<Sync> {
	void lock() {}
	void unlock() {}
};

template<>
struct Mutex<Async> {
	void lock() {
		mtx_.lock();
	}
	void unlock() {
		mtx_.unlock();
	}
private:
	std::mutex mtx_;
};