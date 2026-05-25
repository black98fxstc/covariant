#include <iostream>
#include <string>
#include <fstream>
#include <functional>

// A popular, header-only library for command-line parsing.
// See setup instructions below.
#include <cxxopts.hpp>

// The Covariant class header.
#include "Covariant.hpp"

struct Params {
    std::string filename;
    unsigned dimension;
    unsigned normal;
    unsigned snake;
    unsigned exponential;
    size_t events;
    bool ascii;
};

template <unsigned Dimension>
int do_it(const Params &params) {
    TestData<Dimension> test_data;
    typename TestData<Dimension>::RandomSample test_sample;
    for (unsigned i = 0; i < params.normal; i++)
        test_sample.subpopulation(new typename TestData<Dimension>::Normal());
    for (unsigned i = 0; i < params.snake; i++)
        test_sample.subpopulation(new typename TestData<Dimension>::Snake());
    for (unsigned i = 0; i < params.exponential; i++)
        test_sample.subpopulation(new typename TestData<Dimension>::Exponential());
    test_data.generate(test_sample, params.events);
    std::string ext = params.ascii ? ".txt" : ".dat";
    std::cout << "Saving " << test_data.size() << " events to " << params.filename << ext << "..." << std::endl;
    test_data.write(params.filename + ext, params.ascii);
    return 0;
}

int main(int argc, char *argv[])
{
    Params params;

    // Set up the command-line options parser.
    cxxopts::Options options("WeightyCLI", "Generate test data for the Weighty class");

    // Add the command-line options.
    options.add_options()
        ("f,file", "File name for generated test data", cxxopts::value<std::string>()->default_value("test_data"))
        ("d,dimension", "Dimension of the events", cxxopts::value<unsigned>()->default_value("2"))
        ("n,normal", "Normal distributions", cxxopts::value<unsigned>()->default_value("3"))
        ("exponential", "Exponential distributions", cxxopts::value<unsigned>()->default_value("0"))
        ("snake", "Snake distributions", cxxopts::value<unsigned>()->default_value("0"))
        ("e,events", "Number of events to generate", cxxopts::value<size_t>()->default_value("10000"))
        ("a,ascii", "Use ASCII format for data files", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage");

    // Add the command-line options.
    options.add_options()
        ("f,file", "File name for generated test data", cxxopts::value<std::string>()->default_value("test_data"))
        ("d,dimension", "Dimension of the events", cxxopts::value<unsigned>()->default_value("2"))
        ("n,normal", "Normal distributions", cxxopts::value<unsigned>()->default_value("3"))
        ("exponential", "Exponential distributions", cxxopts::value<unsigned>()->default_value("0"))
        ("snake", "Snake distributions", cxxopts::value<unsigned>()->default_value("0"))
        ("e,events", "Number of events to generate", cxxopts::value<size_t>()->default_value("10000"))
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
    params.normal = result["normal"].as<unsigned>();
    params.snake = result["snake"].as<unsigned>();
    params.exponential = result["exponential"].as<unsigned>();
    params.events = result["events"].as<size_t>();
    params.ascii = result["ascii"].as<bool>();

    std::cout << "Program running with --normal=" << params.normal << " --snake=" << params.snake << " --exponential=" << params.exponential << std::endl;
    std::cout << "Generating " << params.events << " events in " << params.dimension << " dimensions..." << std::endl;

    switch (params.dimension)
    {
    case 2: return do_it<2>(params);
    case 3: return do_it<3>(params);
    case 4: return do_it<4>(params);
    default:
        std::cerr << "Error: Unsupported dimension: " << params.dimension << std::endl;
        return 1;
    }
}