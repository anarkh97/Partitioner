#ifndef __LOG_UTILS_H__
#define __LOG_UTILS_H__

#include <cstdio>

namespace LogUtils
{

double WallTime();

void PrintHeader( );
void Print(const char* format, ...);
void PrintSuccess(const char* format, ...);
void PrintError(const char* format, ...);
void PrintWarning(const char* format, ...);

void PrintNoBarrier(const char* format, ...);
void PrintSuccessNoBarrier(const char* format, ...);
void PrintErrorNoBarrier(const char* format, ...);
void PrintWarningNoBarrier(const char* format, ...);

void Print(FILE*, const char* format, ...);
void PrintSuccess(FILE*, const char* format, ...);
void PrintError(FILE*, const char* format, ...);
void PrintWarning(FILE*, const char* format, ...);

void PrintNoBarrier(FILE*, const char* format, ...);
void PrintSuccessNoBarrier(FILE*, const char* format, ...);
void PrintErrorNoBarrier(FILE*, const char* format, ...);
void PrintWarningNoBarrier(FILE*, const char* format, ...);

};  // namespace LogUtils

#endif
