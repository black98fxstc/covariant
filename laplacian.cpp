#include <iostream>
#include <string>
#include <fstream>
#include <cxxopts.hpp>
#include <functional>

#include "Covariant.hpp"

struct Params
{
    std::string filename;
    unsigned dimension;
    unsigned grid;
    float smooth;
    float threshold;
    bool visual;
    bool verbose;
    bool antialias;
    bool verify;
    bool grow;
    bool ascii;
};

template <unsigned Dimension>
int do_it(const Params &params, unsigned grid_size)
{
    TestData<Dimension> events;
    std::string ext = params.ascii ? ".txt" : ".dat";
    std::cout << "Loading events from " << params.filename << ext << "..." << std::endl;
    if (!events.read(params.filename + ext, params.ascii))
    {
        std::cerr << "Error: Could not open event file for loading: " << params.filename + ext << std::endl;
        return 1;
    }
    std::cout << "Loaded " << events.size() << " events." << std::endl;

    std::cout << "Processing " << events.size() << " events..." << std::endl;
    Laplacian<Dimension> laplace(grid_size, true);
    laplace.visualize = params.visual;
    laplace.antialias = params.antialias;
    laplace.verify = params.verify;
    laplace.verbose = params.verbose;
    for (const auto &e : events)
        laplace.event(e);
    std::cout << "Found " << events.size() << " valid events..." << std::endl;

    std::cout << "Analyzing the sample..." << std::endl;
    laplace.analyze(params.smooth, params.threshold);
    if (laplace.differentialError() > .0001)
        std::cout << "Differential equation solution is unusually bad " << laplace.differentialError() << std::endl;
    else
        std::cout << "Consistency checks passed..." << std::endl;

    std::cout << "Performing Laplacian clustering..." << std::endl;
    unsigned found = laplace.cluster(params.threshold, params.grow);
    std::vector<unsigned short> classes;
    classes.reserve(events.size());
    for (const auto &e : events)
        classes.push_back(laplace.classify(e));
    std::cout << "Found " << found << " clusters." << std::endl;

    if (params.ascii)
    {
        std::ofstream out(params.filename + ".cluster");
        if (!out)
        {
            std::cerr << "Error: Could not open cluster file for writing: " << params.filename + ".cluster" << std::endl;
            return 1;
        }
        for (auto c : classes) out << c << "\n";
        out.close();
    }
    else
    {
        std::ofstream out(params.filename + ".cluster", std::ios::binary | std::ios::trunc);
        if (!out)
        {
            std::cerr << "Error: Could not open cluster file for writing: " << params.filename + ".cluster" << std::endl;
            return 1;
        }
        out.write(reinterpret_cast<const char *>(classes.data()), classes.size() * sizeof(unsigned short));
        out.close();
    }
    std::cout << "Saved to file " << params.filename + ".cluster" << std::endl;

    return 0;
}

int main(int argc, char *argv[])
{
    Params params;

    // Set up the command-line options parser.
    cxxopts::Options options("CovariantCLI", "A command-line interface for the Laplacian clustering algorithm");

    // Add the command-line options.
    options.add_options()
        ("f,file", "File name for data", cxxopts::value<std::string>()->default_value("test_data"))
        ("d,dimension", "Dimension of the events", cxxopts::value<unsigned>()->default_value("2"))
        ("g,grid", "Grid resolution", cxxopts::value<unsigned>())
        ("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("0.01"))
        ("t,threshold", "Consistency threshold", cxxopts::value<float>()->default_value("0.001"))
        ("visual", "Save files for visualization", cxxopts::value<bool>()->default_value("false"))
        ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
        ("antialias", "Enable antialiasing (developer feature)", cxxopts::value<bool>()->default_value("true"))
        ("verify", "Enable consistency verification", cxxopts::value<bool>()->default_value("true"))
        ("grow", "Unclaimed events are assigned to the nearest cluster", cxxopts::value<bool>()->default_value("false"))
        ("a,ascii", "Use ASCII format for data files", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage");

    options.parse_positional({"file"});

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0;
    }

    params.filename = result["file"].as<std::string>();
    params.dimension = result["dimension"].as<unsigned>();
    params.smooth = result["smooth"].as<float>();
    params.threshold = result["threshold"].as<float>();
    params.visual = result["visual"].as<bool>();
    params.antialias = result["antialias"].as<bool>();
    params.verify = result["verify"].as<bool>();
    params.grow = result["grow"].as<bool>();
    params.verbose = result["verbose"].as<bool>();
    params.ascii = result["ascii"].as<bool>();

    if (result.count("grid"))
    {
        params.grid = result["grid"].as<unsigned>();
    } 
    else 
    {
        switch (params.dimension) {
            case 2:  params.grid = 256; break;
            case 3:  params.grid = 64;  break;
            case 4:  params.grid = 32;  break;
            default: params.grid = 32;  break;
        }
    }

    std::cout << "Program running with dimension=" << params.dimension << " grid=" << params.grid 
              << " filename=" << params.filename << " smooth=" << params.smooth << " threshold=" << params.threshold 
              << " antialias=" << (params.antialias ? "on" : "off") 
              << " verify=" << (params.verify ? "on" : "off")
              << " grow=" << (params.grow ? "on" : "off") << std::endl;
    switch (params.dimension)
    {
    case 2:
        return do_it<2>(params, params.grid);
    case 3:
        return do_it<3>(params, params.grid);
    case 4:
        return do_it<4>(params, params.grid);
    default:
        std::cerr << "Error: Unsupported dimension: " << params.dimension << std::endl;
        return 1;
    }
};