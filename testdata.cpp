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
    std::string file_param;
    unsigned dimension_param, normal_param, snake_param, exponential_param;
    size_t events_param;
    bool ascii_param;

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

    options.parse_positional({"file"});

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0;
    }

    file_param = result["file"].as<std::string>();
    dimension_param = result["dimension"].as<unsigned>();
    normal_param = result["normal"].as<unsigned>();
    snake_param = result["snake"].as<unsigned>();
    exponential_param = result["exponential"].as<unsigned>();
    events_param = result["events"].as<size_t>();
    ascii_param = result["ascii"].as<bool>();

    std::cout << "Program running with --normal=" << normal_param << " --snake=" << snake_param << " --exponential=" << exponential_param << std::endl;
    std::cout << "Generating " << events_param << " events in " << dimension_param << " dimensions..." << std::endl;

    switch (dimension_param)
    {
    case 2:
    {
        TestData<2> test_data;
        TestData<2>::RandomSample test_sample;
        for (unsigned i = 0; i < normal_param; i++)
            test_sample.subpopulation(new TestData<2>::Normal());
        for (unsigned i = 0; i < snake_param; i++)
            test_sample.subpopulation(new TestData<2>::Snake());
        for (unsigned i = 0; i < exponential_param; i++)
            test_sample.subpopulation(new TestData<2>::Exponential());
        test_data.generate(test_sample, events_param);
        std::string ext = ascii_param ? ".txt" : ".dat";
        std::cout << "Saving " << test_data.size() << " events to " << file_param << ext << "..." << std::endl;
        test_data.write(file_param + ext, ascii_param);
    }
    break;
    case 3:
    {
        TestData<3> test_data;
        TestData<3>::RandomSample test_sample;
        for (unsigned i = 0; i < normal_param; i++)
            test_sample.subpopulation(new TestData<3>::Normal());
        for (unsigned i = 0; i < snake_param; i++)
            test_sample.subpopulation(new TestData<3>::Snake());
        for (unsigned i = 0; i < exponential_param; i++)
            test_sample.subpopulation(new TestData<3>::Exponential());
        test_data.generate(test_sample, events_param);
        std::string ext = ascii_param ? ".txt" : ".dat";
        std::cout << "Saving " << test_data.size() << " events to " << file_param << ext << "..." << std::endl;
        test_data.write(file_param + ext, ascii_param);
    }
    break;
    default:
        std::cerr << "Error: Unsupported dimension: " << dimension_param << std::endl;
        return 1;
    }
}