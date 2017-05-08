#include <stdint.h>
#ifdef _MSC_VER
#else
#include <time.h>
#endif

struct High_res_timer {
	
	High_res_timer() {
#ifdef _MSC_VER
#else
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time_);
#endif
	}

	uint64_t nanoseconds() {
#ifdef _MSC_VER
		return 0;		
#else
		struct timespec end_time;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time);
		return end_time.tv_nsec - time_.tv_nsec;
#endif
	}

	double microseconds()
	{
		return nanoseconds() / 1000.0;
	}

private:
#ifdef _MSC_VER
#else
	timespec time_;
#endif

};