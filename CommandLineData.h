#ifndef __COMMAND_LINE_DATA_H__
#define __COMMAND_LINE_DATA_H__

#include <mpi.h>
#include <string>
#include <petsc.h>
#include <argparse/argparse.hpp>

#include "LogUtils.h"

enum PatchPartitionBackend { METIS=0, JET };

struct CommandLineData {

    argparse::ArgumentParser cli;

    //!< Input
    std::string mesh_file;

    //!< Patch partitioning
    std::string backend_name = "metis";
    PatchPartitionBackend backend = PatchPartitionBackend::METIS;
    int npatches  = 0;
    int max_bytes_patch = 0;
    int max_elems_patch = 0;
    int seed = 0;

    //!< Output
    std::string prefix;
    std::string out_file;

    explicit CommandLineData(MPI_Comm& comm, int argc, char* argv[]);
    ~CommandLineData() { };
};

//-----------------------------------------------------------------------------

inline
CommandLineData::CommandLineData(MPI_Comm& comm, int argc, char* argv[])
    : cli("npart", "1.0", argparse::default_arguments::none)
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

    //!<---------------------------------------------------
    //!< Help and version
    //!<---------------------------------------------------
    cli.add_argument("-h", "--help")
        .help("shows help message and exits")
        .flag();

    cli.add_argument("-v", "--version")
        .help("prints version information and exits")
        .flag();

    //!<---------------------------------------------------
    //!< Input
    //!<---------------------------------------------------
    cli.add_argument("--mesh-file")
        .help("Path to the input mesh file to partition "
              "(supported types: .msh GMSH v2.2).")
        .required()
        .store_into(mesh_file);

    //!<---------------------------------------------------
    //!< Patch partitioning
    //!<---------------------------------------------------
    cli.add_argument("--backend")
        .help("Patch partitioner (metis or jet).")
        .choices("metis", "jet")
        .default_value("metis")
        .store_into(backend_name);

    //!< Here true means this group is required.
    auto& patch_size = cli.add_mutually_exclusive_group(true);

    patch_size.add_argument("--npatches")
        .help("Number of patches per rank.")
        .store_into(npatches);

    patch_size.add_argument("--bytes-per-patch")
        .help("Per-patch distinct DOF byte cap, this sets the patch count.")
        .store_into(max_bytes_patch);

    cli.add_argument("--elems-per-patch")
        .help("Maximum elements per patch (0 means no bound).")
        .required()
        .default_value(0)
        .store_into(max_elems_patch);

    cli.add_argument("--seed")
        .help("Partitioner random seed.")
        .default_value(0)
        .store_into(seed);

    //!<---------------------------------------------------
    //!< Output
    //!<---------------------------------------------------
    cli.add_argument("--prefix")
        .help("Directory for output files.")
        .default_value("./")
        .store_into(prefix);

    cli.add_argument("--out-file")
        .help("Name of the output file (partition and patch numbers will be "
              "appended to this name).")
        .required()
        .store_into(out_file);

    //!<---------------------------------------------------
    //!< Parse
    //!<---------------------------------------------------
    try {
        cli.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        if (cli.is_used("--help")) {
            LogUtils::Print("%s", cli.help().str().c_str());
            PetscFinalize();
            MPI_Abort(comm,0);
        }
        LogUtils::PrintError("%s\n\n%s\nRun with --help for details.",
                             err.what(), cli.usage().c_str());
    }

    if (cli.is_used("--help")) {
        LogUtils::Print("%s", cli.help().str().c_str());
        PetscFinalize();
        MPI_Abort(comm,0);
    }
    if (cli.is_used("--version")) {
        LogUtils::Print("npart v1.0.0");
        PetscFinalize();
        MPI_Abort(comm,0);
    }

    //!<---------------------------------------------------
    //!< Validate
    //!<---------------------------------------------------
    if (backend_name == "metis") {
        backend = PatchPartitionBackend::METIS;
    }
    else if (backend_name == "jet") {
#ifdef HAVE_JET
        backend = PatchPartitionBackend::JET;
#else
        LogUtils::PrintError("--backend jet requested, "
                             "but npart was built without JET.");
#endif
    }

    if (cli.is_used("--npatches") && npatches <= 0) {
        LogUtils::PrintError("--npatches must be positive (got %d).",
                             npatches);
    }
    if (cli.is_used("--bytes-per-patch") && max_bytes_patch <= 0) {
        LogUtils::PrintError("--bytes-per-patch must be positive (got %d).",
                             max_bytes_patch);
    }
    if (max_elems_patch < 0) {
        LogUtils::PrintError("--elems-per-patch must be non-negative (got %d).",
                             max_elems_patch);
    }

    out_file = prefix + "/" + out_file + ".vtu";
}

#endif
