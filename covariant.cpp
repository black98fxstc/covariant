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
    unsigned dimension_param, cluster_param;
    float smooth_param, threshold_param;
    bool visual_param, verbose_param;

    // Set up the command-line options parser.
    cxxopts::Options options("WeightyCLI", "Generate test data for the Weighty class");

    // Add the command-line options.
    options.add_options()
        ("f,file", "Output filename for generated data", cxxopts::value<std::string>()->default_value("test_data"))
        ("c,cluster", "Process data from the specified cluster", cxxopts::value<unsigned>()->default_value("0"))
        ("d,dimension", "Dimension of the events", cxxopts::value<unsigned>()->default_value("2"))
        ("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("0.01"))
        ("t,threshold", "Consistency threshold", cxxopts::value<float>()->default_value("0.001"))
        ("visual", "Enable visualization", cxxopts::value<bool>()->default_value("false"))
        ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
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
    cluster_param = result["dimension"].as<unsigned>();
    smooth_param = result["smooth"].as<float>();
    threshold_param = result["threshold"].as<float>();
    visual_param = result["visual"].as<bool>();
    verbose_param = result["verbose"].as<bool>();

    std::cout << "Program running with dimension=" << dimension_param << " filename=" << filename_param << "smooth=" << smooth_param << " --threshold=" << threshold_param << std::endl;
    std::cout << "Analyzing events in " << dimension_param << " dimensions..." << std::endl;

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

        Covariant<2> covariant(256, true);
        covariant.visualize = visual_param;
        covariant.verbose = verbose_param;
        size_t valid_events = 0;
        if (result.count("cluster"))
        {
            std::vector<unsigned short> classes(events.size());
            std::ifstream in(filename_param + ".cluster", std::ios::binary);
            in.read(reinterpret_cast<char *>(classes.data()), events.size() * sizeof(unsigned short));
            in.close();

            std::cout << "Processing cluster " << cluster_param << " of " << events.size() << " events..." << std::endl;
            for (unsigned i = 0; i < events.size(); i++)
            {
                if (classes[i] != cluster_param)
                    continue;
                if (covariant.event(events[i]))
                    valid_events++;
            }
        }
        else
        {
            std::cout << "Processing " << events.size() << " events..." << std::endl;
            for (const auto &e : events)
                if (covariant.event(e))
                    valid_events++;
        }
        std::cout << "Found " << events.size() << " events..." << std::endl;

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

        if (!result.count("cluster"))
        {
            std::cout << "Performing Laplacian clustering..." << std::endl;
            unsigned found = covariant.cluster(threshold_param);
            std::vector<unsigned short> classes(events.size());
            for (const auto& e : events)
                classes.push_back(covariant.classify(e));
            std::cout << "Found " << found << " clusters..." << std::endl;

            std::ofstream out(filename_param + ".cluster", std::ios::binary | std::ios::trunc);
            out.write(reinterpret_cast<const char *>(classes.data()), covariant.size() * sizeof(unsigned short));
            out.close();
            std::cout << "Saved to file " << filename_param + ".cluster" << std::endl;
        }
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