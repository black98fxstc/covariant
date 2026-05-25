#include <iostream>
#include <string>
#include <fstream>

// A popular, header-only library for command-line parsing.
// See setup instructions below.
#include <cxxopts.hpp>
#include <functional>

// The Covariant class header.
// #include "Weighty.hpp"
// #include "Weighty.hpp"
// #include "Laplacian.hpp"
#include "Covariant.hpp"

const unsigned Dimension = 2;

int main(int argc, char *argv[])
{
    std::string file_param;
    unsigned dimension_param;
    float smooth_param, threshold_param;
    bool visual_param, verbose_param, ascii_param;

    // Set up the command-line options parser.
    cxxopts::Options options("CovariantCLI", "A command-line interface for the Laplacian clustering algorithm");

    // Add the command-line options.
    options.add_options()
    ("f,file", "File name for data", cxxopts::value<std::string>()->default_value("test_data"))
    ("d,dimension", "Dimension of the events", cxxopts::value<unsigned>()->default_value("2"))
    ("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("0.01"))
    ("t,threshold", "Consistency threshold", cxxopts::value<float>()->default_value("0.001"))
    ("visual", "Enable visualization", cxxopts::value<bool>()->default_value("false"))
    ("v,verbose", "Verbose output", cxxopts::value<bool>()->default_value("false"))
    ("a,ascii", "Use ASCII format for data files", cxxopts::value<bool>()->default_value("false"))
    ("h,help", "Print usage");

    options.parse_positional({"file"});

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0;
    }

    file_param = result["file"].as<std::string>();
    dimension_param = result["dimension"].as<unsigned>();
    smooth_param = result["smooth"].as<float>();
    threshold_param = result["threshold"].as<float>();
    visual_param = result["visual"].as<bool>();
    verbose_param = result["verbose"].as<bool>();
    ascii_param = result["ascii"].as<bool>();

    std::cout << "Program running with dimension=" << dimension_param << " filename=" << file_param << " smooth=" << smooth_param << " threshold=" << threshold_param << std::endl;

    switch (dimension_param)
    {
    case 2:
    {
        TestData<2> events;
        std::string ext = ascii_param ? ".txt" : ".dat";
        std::cout << "Loading events from " << file_param << ext << "..." << std::endl;
        if (!events.read(file_param + ext, ascii_param))
        {
            std::cerr << "Error: Could not open event file for loading: " << file_param + ext << std::endl;
            return 1;
        }
        std::cout << "Loaded " << events.size() << " events." << std::endl;

        std::cout << "Processing " << events.size() << " events..." << std::endl;
        Laplacian<2> laplace(256, true);
        laplace.visualize = visual_param;
        laplace.verbose = verbose_param;
        for (const auto &e : events)
            laplace.event(e);

        std::cout << "Analyzing the sample..." << std::endl;

        laplace.analyze(smooth_param, threshold_param);
        if (laplace.differentialEquation() > .0001)
        {
            std::cout << "Differential equation solution is unusually bad " << laplace.differentialEquation() << std::endl;
        }
        else
        {
            std::cout << "Consistency checkes passed..." << std::endl;
        }

        std::cout << "Performing Laplacian clustering..." << std::endl;
        unsigned found = laplace.cluster(threshold_param);
        std::vector<unsigned short> classes(events.size());
        // for (const auto& e : events)
        //     classes.push_back(laplace.classify(e));
        std::cout << "Found " << found << " clusters." << std::endl;
    }
    break;

    case 3:
    {
        TestData<3> events;
        std::string ext = ascii_param ? ".txt" : ".dat";
        std::cout << "Loading events from " << file_param << ext << "..." << std::endl;
        if (!events.read(file_param + ext, ascii_param))
        {
            std::cerr << "Error: Could not open event file for loading: " << file_param + ext << std::endl;
            return 1;
        }
        std::cout << "Loaded " << events.size() << " events." << std::endl;

        std::cout << "Processing " << events.size() << " events..." << std::endl;
        Laplacian<3> laplace(64, true);
        laplace.visualize = visual_param;
        laplace.verbose = verbose_param;
        for (const auto &e : events)
            laplace.event(e);

        std::cout << "Analyzing the sample..." << std::endl;

        laplace.analyze(smooth_param, threshold_param);
        if (laplace.differentialEquation() > .0001)
        {
            std::cout << "Differential equation solution is unusually bad " << laplace.differentialEquation() << std::endl;
        }
        else
        {
            std::cout << "Consistency checkes passed..." << std::endl;
        }

        std::cout << "Performing Laplacian clustering..." << std::endl;
        unsigned found = laplace.cluster(threshold_param);
        std::vector<unsigned short> classes(events.size());
        // for (const auto& e : events)
        //     classes.push_back(laplace.classify(e));
        std::cout << "Found " << found << " clusters." << std::endl;
    }
    break;

    case 4:
    {
        TestData<4> events;
        std::string ext = ascii_param ? ".txt" : ".dat";
        std::cout << "Loading events from " << file_param << ext << "..." << std::endl;
        if (!events.read(file_param + ext, ascii_param))
        {
            std::cerr << "Error: Could not open event file for loading: " << file_param + ext << std::endl;
            return 1;
        }
        std::cout << "Loaded " << events.size() << " events." << std::endl;

        std::cout << "Processing " << events.size() << " events..." << std::endl;
        Laplacian<4> laplace(32, true);
        laplace.visualize = visual_param;
        laplace.verbose = verbose_param;
        for (const auto &e : events)
            laplace.event(e);

        std::cout << "Analyzing the sample..." << std::endl;

        laplace.analyze(smooth_param, threshold_param);
        if (laplace.differentialEquation() > .0001)
        {
            std::cout << "Differential equation solution is unusually bad " << laplace.differentialEquation() << std::endl;
        }
        else
        {
            std::cout << "Consistency checkes passed..." << std::endl;
        }

        std::cout << "Performing Laplacian clustering..." << std::endl;
        unsigned found = laplace.cluster(threshold_param);
        std::vector<unsigned short> classes(events.size());
        // for (const auto& e : events)
        //     classes.push_back(laplace.classify(e));
        std::cout << "Found " << found << " clusters." << std::endl;
    }
    break;

    break;
    default:
        std::cerr << "Error: Unsupported dimension: " << dimension_param << std::endl;
        return 1;
    }
};