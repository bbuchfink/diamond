#include <stdint.h>
#ifdef _MSC_VER
#else
#include <time.h>
#endif

struct High_res_timer {
	
	High_res_timer():
#ifdef _MSC_VER
	time_(__rdtsc())
#else
	time_(0)
#endif
	{
	}

	uint64_t get() const
	{
#ifdef _MSC_VER
		return __rdtsc() - time_;
#else
		return 0;
#endif
	}

	uint64_t nanoseconds() {
#ifdef _MSC_VER
		return 0;		
#else
		return 0;
#endif
	}

	double microseconds()
	{
		return nanoseconds() / 1000.0;
	}

private:
	unsigned long long time_;

};