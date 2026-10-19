#include <LogUtils.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <cstdlib>

extern MPI_Comm comm;

namespace LogUtils
{

double WallTime()
{
    return MPI_Wtime();
}

const char* RESET = "\033[0m";
const char* EMPTY = "";
const char* RED = "\033[31m*** Error: ";
const char* YELLOW = "\033[33m*** Warning: ";
const char* GREEN = "\033[32m";

void VPrint(FILE* file, const char* prefix, const char* suffix,
            const char* format, va_list& argsp)
{
    fprintf(file, "%s", prefix);
    vfprintf(file, format, argsp);
    fprintf(file, "%s\n", suffix);
}

void PrintHeader( )
{
    // custom header here.
}

void Print(const char* format, ...)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
        va_list Argp;
        va_start(Argp, format);

        VPrint(stdout, EMPTY, EMPTY, format, Argp);

        va_end(Argp);
    }
    MPI_Barrier(comm);
    return;
}

void PrintSuccess(const char* format, ...)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
        va_list Argp;
        va_start(Argp, format);

        VPrint(stdout, GREEN, RESET, format, Argp);

        va_end(Argp);
    }
    MPI_Barrier(comm);
    return;
}

void PrintError(const char* format, ...)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
        va_list Argp;
        va_start(Argp, format);

        VPrint(stderr, RED, RESET, format, Argp);

        va_end(Argp);
    }
    MPI_Abort(comm, -1);
    return;
}

void PrintWarning(const char* format, ...)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0)
    {
        va_list Argp;
        va_start(Argp, format);

        VPrint(stdout, YELLOW, RESET, format, Argp);

        va_end(Argp);
    }
    MPI_Barrier(comm);
    return;
}

void PrintNoBarrier(const char* format, ...)
{
    //int rank;
    //MPI_Comm_rank(comm, &rank);
    //if (!rank)
    //{
        va_list Argp;
        va_start(Argp, format);

        VPrint(stdout, EMPTY, EMPTY, format, Argp);

        va_end(Argp);
    //}
    //MPI_Barrier(comm);
    return;
}

void PrintSuccessNoBarrier(const char* format, ...)
{
    //int rank;
    //MPI_Comm_rank(comm, &rank);
    //if (!rank)
    //{
        va_list Argp;
        va_start(Argp, format);

        VPrint(stdout, GREEN, RESET, format, Argp);

        va_end(Argp);
    //}
    //MPI_Barrier(comm);
    return;
}

void PrintErrorNoBarrier(const char* format, ...)
{
    //int rank;
    //MPI_Comm_rank(comm, &rank);
    //if (!rank)
    //{
        va_list Argp;
        va_start(Argp, format);

        VPrint(stderr, RED, RESET, format, Argp);

        va_end(Argp);
    //}
    MPI_Abort(comm, -1);
}

void PrintWarningNoBarrier(const char* format, ...)
{
    //int rank;
    //MPI_Comm_rank(comm, &rank);
    //if (!rank)
    //{
        va_list Argp;
        va_start(Argp, format);

        VPrint(stdout, YELLOW, RESET, format, Argp);

        va_end(Argp);
    //}
    //MPI_Barrier(comm);
    return;
}


void Print(FILE* file, const char* format, ...)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
        va_list Argp;
        va_start(Argp, format);

        VPrint(file, EMPTY, EMPTY, format, Argp);

        va_end(Argp);
    }
    MPI_Barrier(comm);
    return;
}

void PrintSuccess(FILE* file, const char* format, ...)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
        va_list Argp;
        va_start(Argp, format);

        VPrint(file, GREEN, RESET, format, Argp);

        va_end(Argp);
    }
    MPI_Barrier(comm);
    return;
}

void PrintError(FILE* file, const char* format, ...)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0) {
        va_list Argp;
        va_start(Argp, format);

        VPrint(file, RED, RESET, format, Argp);

        va_end(Argp);
    }
    MPI_Abort(comm, -1);
    return;
}

void PrintWarning(FILE* file, const char* format, ...)
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    if (rank == 0)
    {
        va_list Argp;
        va_start(Argp, format);

        VPrint(file, YELLOW, RESET, format, Argp);

        va_end(Argp);
    }
    MPI_Barrier(comm);
    return;
}

void PrintNoBarrier(FILE* file, const char* format, ...)
{
    //int rank;
    //MPI_Comm_rank(comm, &rank);
    //if (!rank)
    //{
        va_list Argp;
        va_start(Argp, format);

        VPrint(file, EMPTY, EMPTY, format, Argp);

        va_end(Argp);
    //}
    //MPI_Barrier(comm);
    return;
}

void PrintSuccessNoBarrier(FILE* file, const char* format, ...)
{
    //int rank;
    //MPI_Comm_rank(comm, &rank);
    //if (!rank)
    //{
        va_list Argp;
        va_start(Argp, format);

        VPrint(file, GREEN, RESET, format, Argp);

        va_end(Argp);
    //}
    //MPI_Barrier(comm);
    return;
}

void PrintErrorNoBarrier(FILE* file, const char* format, ...)
{
    //int rank;
    //MPI_Comm_rank(comm, &rank);
    //if (!rank)
    //{
        va_list Argp;
        va_start(Argp, format);

        VPrint(file, RED, RESET, format, Argp);

        va_end(Argp);
    //}
    MPI_Abort(comm, -1);
}

void PrintWarningNoBarrier(FILE* file, const char* format, ...)
{
    //int rank;
    //MPI_Comm_rank(comm, &rank);
    //if (!rank)
    //{
        va_list Argp;
        va_start(Argp, format);

        VPrint(file, YELLOW, RESET, format, Argp);

        va_end(Argp);
    //}
    //MPI_Barrier(comm);
    return;
}

}  // namespace LogUtils
