#ifndef _FILESTACK_H_
#define _FILESTACK_H_

#include <string>
#ifndef WIN32
#include <unistd.h>
#include <fcntl.h>
#endif


#ifdef DEBUG
#define DBG(x) (std::cerr << __PRETTY_FUNCTION__ << ":" << __LINE__ << " " << x << std::endl)
#else
#define DBG(x)
#endif


class FileStack {
    public:
        FileStack();
        FileStack(const std::string & file_name);
        FileStack(const std::string & file_name, int maximum_line_length);
        ~FileStack();

        size_t size();

        int pop(int & i);
        int pop(std::string & buf);
        int pop(std::string & buf, size_t & size_after_pop);
        int pop_non_locked(std::string & buf);

        int top(int & i);
        int top(std::string & buf);

        int remove(const std::string & line);

        int64_t push(int i);
        int64_t push(const std::string & buf);
        int64_t push(const std::string & buf, size_t & size_after_push);
        int64_t push_non_locked(const std::string & buf);

        int get_max_line_length();
        int set_max_line_length(int n);

        int clear();

        int lock();
        int unlock();

        bool poll_query(const std::string & query, const double sleep_s=0.5, const size_t max_iter=7200);
        bool poll_size(const size_t size, const double sleep_s=0.5, const size_t max_iter=7200);

    private:
        int fd;
        bool locked;
#ifndef WIN32
        struct flock lck;
#endif
        std::string file_name;
        off_t max_line_length;

        int pop(std::string &, const bool, size_t &);
        int pop_non_locked(std::string &, const bool, size_t &);
};

#endif
