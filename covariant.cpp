#include <iostream>
#include <string>
#include <fstream>
#include <functional>

// A popular, header-only library for command-line parsing.
// See setup instructions below.
#include <cxxopts.hpp>

// The Covariant class header.
#include "Covariant.hpp"

int main(int argc, char *argv[])
{
    std::string filename_param;
    unsigned dimension_param;
    float smooth_param, threshold_param;

    // Set up the command-line options parser.
    cxxopts::Options options("WeightyCLI", "Generate test data for the Weighty class");

    // Add the command-line options.
    options.add_options()
        ("f,file", "Output filename for generated data", cxxopts::value<std::string>()->default_value("test_data"))
        ("d,dimension", "Dimension of the events", cxxopts::value<unsigned>()->default_value("2"))
        ("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("0.01"))
        ("t,threshold", "Consistency threshold", cxxopts::value<float>()->default_value("0.001"))
        ("h,help", "Print usage");

    options.parse_positional({"file"});

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0;
    }

    filename_param = result["file"].as<std::string>();
    dimension_param = result["dimension"].as<unsigned>();
    smooth_param = result["smooth"].as<float>();
    threshold_param = result["threshold"].as<float>();

    std::cout << "Program running with dimension=" << dimension_param << " filename=" << filename_param << "smooth=" << smooth_param << " --threshold=" << threshold_param << std::endl;
    std::cout << "Generating events in " << dimension_param << " dimensions..." << std::endl;

    switch (dimension_param)
    {
    case 2:
    {
        TestData<2> events;
        std::cout << "Loading events from " << filename_param << "..." << std::endl;
        // std::ifstream infile(filename_param + ".dat", std::ios::binary);
        if (!events.load(filename_param + ".dat"))
        {
            std::cerr << "Error: Could not open event file for loading: " << filename_param + ".dat" << std::endl;
            return 1;
        }
        std::cout << "Loaded " << events.size() << " events." << std::endl;

        std::cout << "Processing " << events.size() << " events..." << std::endl;
        Covariant<2> covariant(256, true);
        for (const auto &e : events)
            covariant.event(e);

        std::cout << "Analyzing the sample..." << std::endl;

        covariant.analyse(smooth_param, threshold_param);
        if (covariant.factorProbability() > .0001)
        {
            std::cout << "Probability factoring is unusually bad " << covariant.factorProbability() << std::endl;
        }
        else if (covariant.differentialEquation() > .0001)
        {
            std::cout << "Differential equation solution is unusually bad " << covariant.differentialEquation() << std::endl;
        }
        else
        {
            std::cout << "Consistency checkes passed..." << std::endl;
        }

        std::cout << "Performing Laplacian clustering..." << std::endl;
        unsigned found = covariant.cluster(threshold_param);
        std::vector<unsigned short> classes(events.size());
        // for (const auto& e : events)
        //     classes.push_back(laplace.classify(e));
        std::cout << "Found " << found << " clusters." << std::endl;
    }
    break;
    case 3:
    {
    }
    break;
    default:
        std::cerr << "Error: Unsupported dimension: " << dimension_param << std::endl;
        return 1;
    }
};