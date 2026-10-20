#ifndef __COMMAND_LINE_DATA_H__
#define __COMMAND_LINE_DATA_H__

#include <mpi.h>
#include <argparse/argparse.hpp>
#include <string>

struct CommandLineData {

    argparse::ArgumentParser cli;

    int npatches;
    
    std::string mesh_file;

    explicit CommandLineData(MPI_Comm& comm, int argc, char* argv[]);
    ~CommandLineData() { };
};

inline
CommandLineData::CommandLineData(MPI_Comm& comm, int argc, char* argv[])
{

    cli.add_description(
        "                   /$$$$$$$                       /$$    \n"  
        "                  | $$__  $$                     | $$    \n"
        " /$$$$$$$         | $$  \\ $$ /$$$$$$   /$$$$$$  /$$$$$$  \n"
        "| $$__  $$ /$$$$$$| $$$$$$$/|____  $$ /$$__  $$|_  $$_/  \n"
        "| $$  \\ $$|______/| $$____/  /$$$$$$$| $$  \\__/  | $$    \n"
        "| $$  | $$        | $$      /$$__  $$| $$        | $$ /$$\n"
        "| $$  | $$        | $$     |  $$$$$$$| $$        |  $$$$/\n"
        "|__/  |__/        |__/      \\_______/|__/         \\___/  \n"
        "\n"
        "A simple partitioning program for testing interface splits.\n"
    );

    cli.add_argument("mesh_file")
        .help("Path to the input mesh file to parition.")
        .required();

    cli.add_argument("--npatches")
        .help("The number of patches per rank.")
        .required();

    try {
        cli.parse_args(argc, argv);
    }
    catch (const std::exception& err) {

        int rank;
        MPI_Comm_rank(comm, &rank);

        if (rank == 0) {
            fprintf(stdout, "%s", cli.help().str().c_str());
        }

        MPI_Abort(comm,0);
        
    }

    MPI_Barrier(comm);
}

#endif
