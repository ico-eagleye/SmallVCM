#pragma once
#include <cstdarg>
#include <stdio.h>

#define DEBUG_OUTPUT 1
#define DEBUG_X 0
#define DEBUG_Y 0
#define DEBUG_RES_X 512
#define DEBUG_IDX DEBUG_X + DEBUG_Y * DEBUG_RES_X
#define DEBUG_THREAD_ID 0

#if DEBUG_OUTPUT

void dbgPrintf(const char * format, ... ) 
{ 
    va_list args;
    va_start(args, format);
    vfprintf(stdout, format, args);
    fflush(stdout);
};

#define DBG_PRINTFI(idx, format, ...) \
    if (idx && (*idx == DEBUG_IDX) && (DEBUG_THREAD_ID == omp_get_thread_num())) \
    { \
        dbgPrintf(format, __VA_ARGS__); \
    }
#else
inline void dbgPrintf(const char * format, ... ) { };
#define DBG_PRINTFI(idx, format, ...)
#endif





