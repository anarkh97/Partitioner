#ifndef __COMMAND_LINE_DATA_H__
#define __COMMAND_LINE_DATA_H__

#include <mpi.h>
#include <petsc.h>
#include <argparse/argparse.hpp>
#include <exception>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <string>

#include "GlobalTypeDefs.h"
#include "LogUtils.h"

enum PatchPartitionBackend
{
    METIS = 0,
    JET
};

struct CommandLineData
{
    using json = nlohmann::json;

    argparse::ArgumentParser cli;

    //!< Config file
    std::filesystem::path config_path =
        std::filesystem::current_path( ) / "config.json";

    //!< Input
    std::string mesh_file;

    //!< Patch partitioning
    std::string backend_name = "metis";
    PatchPartitionBackend backend = PatchPartitionBackend::METIS;

    int_t npatches = 0;
    int_t max_bytes_patch = 0;
    int_t sm_count = 0;
    int_t shared_bytes_per_sm = 0;
    int_t resident_per_sm = 4;
    int_t min_patches_per_sm = 2;
    int_t reserve_bytes = 0;
    int_t max_elems_patch = 0;
    int_t seed = 0;

    //!< Output
    std::string prefix;
    std::string out_file;

    explicit CommandLineData(MPI_Comm& comm, int argc, char* argv[]);
    ~CommandLineData( ) {};
};

//-----------------------------------------------------------------------------

inline CommandLineData::CommandLineData(MPI_Comm& comm, int argc, char* argv[])
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
        "A simple partitioning program for testing interface splits.\n");

    //!<---------------------------------------------------
    //!< Help and version
    //!<---------------------------------------------------
    cli.add_argument("-h", "--help")
        .help("shows help message and exits")
        .flag( );

    cli.add_argument("-v", "--version")
        .help("prints version information and exits")
        .flag( );

    //!<---------------------------------------------------
    //!< Input
    //!<---------------------------------------------------
    cli.add_argument("--mesh-file")
        .help(
            "Path to the input mesh file to partition "
            "(supported types: .msh GMSH v2.2).")
        .required( )
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

    patch_size.add_argument("--shared-bytes-per-sm")
        .help(
            "Shared memory per SM in bytes, the cap is derived from it "
            "using --resident-per-sm and --reserve-bytes.")
        .store_into(shared_bytes_per_sm);

    cli.add_argument("--sm-count")
        .help(
            "Number of Streaming Multiprocessors (SM) on the target GPU "
            "(0 skips the GPU-fill bound).")
        .default_value(0)
        .store_into(sm_count);

    cli.add_argument("--resident-per-sm")
        .help("Target number of patches resident on one SM at a time.")
        .default_value(4)
        .store_into(resident_per_sm);

    cli.add_argument("--min-patches-per-sm")
        .help("Lower bound on launched patches per SM, to keep every SM busy.")
        .default_value(8)
        .store_into(min_patches_per_sm);

    cli.add_argument("--reserve-bytes")
        .help(
            "Shared memory per patch reserved for tabulations and other "
            "scratch, subtracted from the derived cap.")
        .default_value(0)
        .store_into(reserve_bytes);


    cli.add_argument("--elems-per-patch")
        .help("Maximum elements per patch (0 means no bound).")
        .required( )
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
        .help(
            "Name of the output file (partition and patch numbers will be "
            "appended to this name).")
        .required( )
        .store_into(out_file);

    //!<---------------------------------------------------
    //!< Check for config
    //!<---------------------------------------------------
    if (std::filesystem::exists(config_path))
    {
        LogUtils::Print("Found configuration file %s.",
                        config_path.string( ).c_str( ));

        std::vector<std::string> json_args = {argv[0]};

        std::ifstream file(config_path);

        try
        {
            json conf_data;
            file >> conf_data;

            //!< Input
            if (conf_data.contains("mesh-file"))
            {
                json_args.push_back("--mesh-file");
                json_args.push_back(conf_data["mesh-file"].get<std::string>( ));
            }

            //!< Patch partitioning
            if (conf_data.contains("backend"))
            {
                json_args.push_back("--backend");
                json_args.push_back(conf_data["backend"].get<std::string>( ));
            }
            if (conf_data.contains("npatches"))
            {
                json_args.push_back("--npatches");
                json_args.push_back(
                    std::to_string(conf_data["npatches"].get<int_t>( )));
            }
            if (conf_data.contains("bytes-per-patch"))
            {
                json_args.push_back("--bytes-per-patch");
                json_args.push_back(
                    std::to_string(conf_data["bytes-per-patch"].get<int_t>( )));
            }
            if (conf_data.contains("shared-bytes-per-sm"))
            {
                json_args.push_back("--shared-bytes-per-sm");
                json_args.push_back(std::to_string(
                    conf_data["shared-bytes-per-sm"].get<int_t>( )));
            }
            if (conf_data.contains("sm-count"))
            {
                json_args.push_back("--sm-count");
                json_args.push_back(
                    std::to_string(conf_data["sm-count"].get<int_t>( )));
            }
            if (conf_data.contains("resident-per-sm"))
            {
                json_args.push_back("--resident-per-sm");
                json_args.push_back(
                    std::to_string(conf_data["resident-per-sm"].get<int_t>( )));
            }
            if (conf_data.contains("min-patches-per-sm"))
            {
                json_args.push_back("--min-patches-per-sm");
                json_args.push_back(std::to_string(
                    conf_data["min-patches-per-sm"].get<int_t>( )));
            }
            if (conf_data.contains("reserve-bytes"))
            {
                json_args.push_back("--reserve-bytes");
                json_args.push_back(
                    std::to_string(conf_data["reserve-bytes"].get<int_t>( )));
            }
            if (conf_data.contains("elems-per-patch"))
            {
                json_args.push_back("--elems-per-patch");
                json_args.push_back(
                    std::to_string(conf_data["elems-per-patch"].get<int_t>( )));
            }
            if (conf_data.contains("seed"))
            {
                json_args.push_back("--seed");
                json_args.push_back(
                    std::to_string(conf_data["seed"].get<int_t>( )));
            }

            //!< Output
            if (conf_data.contains("prefix"))
            {
                json_args.push_back("--prefix");
                json_args.push_back(conf_data["prefix"].get<std::string>( ));
            }
            if (conf_data.contains("out-file"))
            {
                json_args.push_back("--out-file");
                json_args.push_back(conf_data["out-file"].get<std::string>( ));
            }

            cli.parse_args(json_args);
        }
        catch (const std::exception& err)
        {
            LogUtils::PrintError("Could not parse config.json. "
                                 "%s\n\n%s\nRun with --help for details.",
                                 err.what( ), cli.usage( ).c_str( ));
        }
    }
    else
    {
        try
        {
            cli.parse_args(argc, argv);
        }
        catch (const std::exception& err)
        {
            if (cli.is_used("--help"))
            {
                LogUtils::Print("%s", cli.help( ).str( ).c_str( ));
                MPI_Abort(comm, 0);
            }
            LogUtils::PrintError("%s\n\n%s\nRun with --help for details.",
                                 err.what( ), cli.usage( ).c_str( ));
        }
    }

    if (cli.is_used("--help"))
    {
        LogUtils::Print("%s", cli.help( ).str( ).c_str( ));
        PetscFinalize( );
        MPI_Abort(comm, 0);
    }
    if (cli.is_used("--version"))
    {
        LogUtils::Print("npart v1.0.0");
        PetscFinalize( );
        MPI_Abort(comm, 0);
    }

    //!<---------------------------------------------------
    //!< Validate
    //!<---------------------------------------------------
    if (backend_name == "metis")
    {
        backend = PatchPartitionBackend::METIS;
    }
    else if (backend_name == "jet")
    {
#ifdef HAVE_JET
        backend = PatchPartitionBackend::JET;
#else
        LogUtils::PrintError(
            "--backend jet requested, "
            "but npart was built without JET.");
#endif
    }

    if (cli.is_used("--npatches") && npatches <= 0)
    {
        LogUtils::PrintError("--npatches must be positive (got %d).", npatches);
    }
    if (cli.is_used("--bytes-per-patch") && max_bytes_patch <= 0)
    {
        LogUtils::PrintError("--bytes-per-patch must be positive (got %d).",
                             max_bytes_patch);
    }
    if (cli.is_used("--shared-bytes-per-sm") && shared_bytes_per_sm < 0)
    {
        LogUtils::PrintError(
            "--shared-bytes-per-sm must be non-negative (got %d).",
            shared_bytes_per_sm);
    }

    if (sm_count < 0)
    {
        LogUtils::PrintError("--sm-count must be non-negative (got %d).",
                             sm_count);
    }
    if (resident_per_sm <= 0)
    {
        LogUtils::PrintError("--resident-per-sm must be positive (got %d).",
                             resident_per_sm);
    }
    if (min_patches_per_sm < resident_per_sm)
    {
        LogUtils::PrintError(
            "--min-patches-per-sm (%d) must be at least "
            "--resident-per-sm (%d).",
            min_patches_per_sm, resident_per_sm);
    }
    if (reserve_bytes < 0)
    {
        LogUtils::PrintError("--reserve-bytes must be non-negative (got %d).",
                             reserve_bytes);
    }
    if (cli.is_used("--shared-bytes-per-sm") &&
        shared_bytes_per_sm / resident_per_sm <= reserve_bytes)
    {
        LogUtils::PrintError("Derived cap is not positive: %d / %d - %d.",
                             shared_bytes_per_sm, resident_per_sm,
                             reserve_bytes);
    }
    if (max_elems_patch < 0)
    {
        LogUtils::PrintError("--elems-per-patch must be non-negative (got %d).",
                             max_elems_patch);
    }

    out_file = prefix + "/" + out_file + ".vtu";
}

#endif
