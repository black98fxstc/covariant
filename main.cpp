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

    // Add the "--normal" option as requested.
    options.add_options()
        ("n,normal", "Normal distributions", cxxopts::value<unsigned>()->default_value("3"));

    // Add a standard help option.
    options.add_options()
        ("h,help", "Print usage");

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    // Retrieve the value of the 'normal' parameter.
    unsigned normal_param = result["normal"].as<unsigned>();
    std::cout << "Program running with --normal=" << normal_param << std::endl;

    // Instantiate Covariant with Dimension=3 and Grid=256 as requested.
    std::cout << "Instantiating Covariant<3, 256>..." << std::endl;
    
    Covariant<3, 256> covariant;

    while (covariant.events() < 1000) {
        Covariant<3, 256>::Event e;
        for (unsigned i = 0; i < covariant.dimension; ++i) {
            e.push_back(static_cast<float>(rand()) / static_cast<float>(RAND_MAX));
        }
        covariant.event(e);
    }
    
    std::cout << "Covariant<3, 256> instantiated successfully." << std::endl;
    std::cout << "Initialized stride values: ";
    for(size_t i = 0; i < covariant.dimension; ++i) {
        std::cout << covariant.stride[i] << (i == covariant.dimension - 1 ? "" : ", ");
    }
    std::cout << std::endl;

    return 0;
}