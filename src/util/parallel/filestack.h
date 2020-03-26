#ifndef _FILESTACK_H_
#define _FILESTACK_H_

#include <iostream>
#include <string>
#include <unistd.h>
#include <fcntl.h>


// #undef DEBUG
#define DEBUG
#ifdef DEBUG
#define DBG(x) (std::cerr << __PRETTY_FUNCTION__ << ":" << __FUNCTION__ << ":" << __LINE__ << " " << x << std::endl)
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

        int pop(std::string & buf);
        int pop(int & i);

        int top(std::string & buf);
        int top(int & i);

        int push(std::string buf);
        int push(int i);

        int get_max_line_length();
        int set_max_line_length(int n);

        int clear();

        int lock();
        int unlock();

        bool poll_query(const std::string & query, const double sleep_s=0.25, const size_t max_iter=1200);
        bool poll_size(const size_t size, const double sleep_s=0.25, const size_t max_iter=1200);

    private:
        int fd;
        bool locked;
        struct flock lck;
        std::string file_name;
        off_t max_line_length;

        int pop(std::string &, const bool);
};

#endif
