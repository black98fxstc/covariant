#include <iostream>
#include <string>

// A popular, header-only library for command-line parsing.
// See setup instructions below.
#include <cxxopts.hpp>

// The Covariant class header.
#include "/home/wmoore/git/covariant/Covariant.hpp"

int main(int argc, char* argv[]) {
    // Set up the command-line options parser.
    cxxopts::Options options("CovariantCLI", "A command-line interface for the Covariant class");

    // Add the command-line options.
    options.add_options()
        ("n,normal", "Normal distributions", cxxopts::value<unsigned>()->default_value("3"))
        ("e,events", "Number of events to generate", cxxopts::value<size_t>()->default_value("10000"))
        ("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("1.0"))
        ("h,help", "Print usage");

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    // Retrieve the values of the parameters.
    unsigned normal_param = result["normal"].as<unsigned>();
    size_t events_param = result["events"].as<size_t>();
    float smooth_param = result["smooth"].as<float>();
    std::cout << "Program running with --normal=" << normal_param << " and --events=" << events_param << " and --smooth=" << smooth_param << std::endl;
    
    Covariant<3> covariant;
    Covariant<3>::Event e;

    while (covariant.events() < events_param) {
        e.clear();
        for (unsigned i = 0; i < covariant.dimension; ++i) {
            e.push_back(static_cast<float>(rand()) / static_cast<float>(RAND_MAX));
        }
        covariant.event(e);
    }
    covariant.parameters(smooth_param);
    
    return 0;
}