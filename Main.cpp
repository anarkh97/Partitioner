#include <mpi.h>
#include <ctime>
#include <cerrno>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <argparse/argparse.hpp>

#include "LogUtils.h"

MPI_Comm comm;
double start_time;

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    start_time = LogUtils::WallTime();
    comm = MPI_COMM_WORLD;




    argparse::ArgumentParser program("test");
    
    program.add_argument("--foo")
        .implicit_value(true)
        .default_value(false);
    
    program.add_argument("--bar");
    
    
    try {
      program.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        LogUtils::PrintError("%s", err.what());
        MPI_Abort(comm,0);
    }
    
    if (program.is_used("--foo")) {
        LogUtils::Print("--foo: %d.", program.get<bool>("--foo") ? 1 : 0); 
    }
    
    if (program.is_used("--bar")) {
        LogUtils::Print("--bar: %s", program.get("--bar").c_str());
    }  

    MPI_Finalize();
    return 0;

}
