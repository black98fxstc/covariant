#include <iostream>
#include <string>
#include <fstream>

// A popular, header-only library for command-line parsing.
// See setup instructions below.
#include <cxxopts.hpp>
#include <functional>

// The Covariant class header.
#include "Covariant.hpp"
#include "TestData.hpp"

typedef float Float;
const unsigned Dimension = 2;

int main(int argc, char* argv[]) {
    // Set up the command-line options parser.
    cxxopts::Options options("CovariantCLI", "A command-line interface for the Covariant class");

    // Add the command-line options.
    options.add_options()
        ("n,normal", "Normal distributions", cxxopts::value<unsigned>()->default_value("3"))
        ("exponential", "Exponential distributions", cxxopts::value<unsigned>()->default_value("0"))
        ("snake", "Snake distributions", cxxopts::value<unsigned>()->default_value("0"))
        ("e,events", "Number of events to generate", cxxopts::value<size_t>()->default_value("10000"))
        ("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("1.0"))
        ("h,help", "Print usage")
        ("save", "Save generated data to a file.", cxxopts::value<std::string>()->implicit_value("covariant.dat"))
        ("load", "Load data from a file instead of generating.", cxxopts::value<std::string>()->implicit_value("covariant.dat"));

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    unsigned normal_param, snake_param, exponential_param;
    size_t events_param;
    float smooth_param;

    normal_param = result["normal"].as<unsigned>();
    snake_param = result["snake"].as<unsigned>();
    exponential_param = result["exponential"].as<unsigned>();
    events_param = result["events"].as<size_t>();
    smooth_param = result["smooth"].as<float>();

    
    TestData<Dimension, Float> events;
    if (result.count("load")) {
        std::string filename = result["load"].as<std::string>();
        std::cout << "Program running with --load=" << filename << " --smooth=" << smooth_param << std::endl;
        std::cout << "Loading events from " << filename << "..." << std::endl;
        std::ifstream infile(filename, std::ios::binary);
        if (!events.load(filename)) {
            std::cerr << "Error: Could not open event file for loading: " << filename << std::endl;
            return 1;
        }
        std::cout << "Loaded " << events.size() << " events." << std::endl;
    } else {
        std::cout << "Program running with --normal=" << normal_param << " --snake=" << snake_param << " --events=" << events_param << " --smooth=" << smooth_param << std::endl;
        std::cout << "Generating " << events_param << " events..." << std::endl;
        TestData<Dimension, Float>::RandomSample test_sample;
        for (unsigned i = 0; i < normal_param; i++)
            test_sample.subpopulation(new TestData<Dimension, Float>::Normal());
        for (unsigned i = 0; i < snake_param; i++)
            test_sample.subpopulation(new TestData<Dimension, Float>::Snake());
        for (unsigned i = 0; i < exponential_param; i++)
            test_sample.subpopulation(new TestData<Dimension, Float>::Exponential());
        events.generate(test_sample, events_param);

        if (result.count("save")) {
            std::string filename = result["save"].as<std::string>();
            std::cout << "Saving " << events.size() << " events to " << filename << "..." << std::endl;
            std::ofstream outfile(filename, std::ios::binary | std::ios::trunc);
            events.save(filename);
        }
    }

    std::cout << "Processing " << events.size() << " events..." << std::endl;
    Covariant<Dimension, Float> covariant(256);
    for (const auto& e : events)
        covariant.event(e);

    std::cout << "Analyzing the sample..." << std::endl;

    covariant.parameters(smooth_param);
    if(covariant.factorProbability() > 0.0001) {
        std::cout << "Probability factoring is unusually bad " << covariant.factorProbability() << std::endl;
    }
    else if(covariant.differentialEquation() > 0.0001) {
        std::cout << "Differential equation solution is unusually bad " << covariant.differentialEquation() << std::endl;
    }
    else {
        std::cout << "Consistency checkes passed..." << std::endl;
    }
    
    std::cout << "Writing data to files..." << std::endl;
    auto write_joint = [&](std::string filename, const Float * data) {
        std::ofstream outfile(filename, std::ios::binary | std::ios::trunc);
        outfile.write(reinterpret_cast<const char*>(&covariant.points), sizeof(covariant.points));
        outfile.write(reinterpret_cast<const char*>(data), covariant.size() * sizeof(Float));
        outfile.close();
    };

    write_joint("w.bin", covariant.w());
    write_joint("f.bin", covariant.f());
    write_joint("QC.bin", covariant.QC());
    write_joint("f1.bin", covariant.f(0));
    write_joint("f2.bin", covariant.f(1));
    write_joint("t11.bin", covariant.t(0, 0));
    write_joint("t12.bin", covariant.t(0, 1));
    write_joint("t21.bin", covariant.t(1, 0));
    write_joint("t22.bin", covariant.t(1, 1));
    write_joint("S1.bin", covariant.S(0));
    write_joint("S2.bin", covariant.S(1));
    write_joint("T1.bin", covariant.T(0));
    write_joint("T2.bin", covariant.T(1));
    write_joint("R.bin", covariant.R());

    auto write_marginal = [&](std::string filename, const int marginal, const Float * data) {
        std::ofstream outfile(filename, std::ios::binary | std::ios::trunc);
        int  m = marginal;
        outfile.write(reinterpret_cast<const char*>(&covariant.points[m]), sizeof(covariant.points[m]));
        outfile.write(reinterpret_cast<const char*>(data), covariant.points[m] * sizeof(Float));
        outfile.close();
    };

    write_marginal("P1.bin", 0, covariant.P(0));
    write_marginal("P2.bin", 1, covariant.P(1));
    write_marginal("Q1.bin", 0, covariant.Q(0));
    write_marginal("Q2.bin", 1, covariant.Q(1));

    std::cout << "Data written successfully." << std::endl;
 
    std::cout << "total R: " << covariant.tot_R() << "  var R: " << covariant.var_R() << "  Significance: " << covariant.tot_R() / std::sqrt(covariant.var_R()) << std::endl;
   
    return 0;
}