#ifndef SYSTEM_C_H_
#define SYSTEM_C_H_

#ifdef _MSC_VER
typedef __int64 ssize_t;
#define STDIN_FILENO _fileno(stdin)
#define STDOUT_FILENO _fileno(stdout)
#endif

#endif