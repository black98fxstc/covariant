#include <iostream>
#include <string>
#include <random>
#include <fstream>

// A popular, header-only library for command-line parsing.
// See setup instructions below.
#include <cxxopts.hpp>
#include <functional>

// The Covariant class header.
#include "Covariant.hpp"

std::mt19937& get_rng() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return gen;
}
static auto& rng = get_rng();

float random_float() {
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    return dist(get_rng());
}

float random_normal(float mean, float stddev) {
    std::normal_distribution<float> dist(mean, stddev);
    return dist(get_rng());
}

int main(int argc, char* argv[]) {
    // Set up the command-line options parser.
    cxxopts::Options options("CovariantCLI", "A command-line interface for the Covariant class");

    // Add the command-line options.
    options.add_options()
        ("n,normal", "Normal distributions", cxxopts::value<unsigned>()->default_value("3"))
        ("e,events", "Number of events to generate", cxxopts::value<size_t>()->default_value("10000"))
        ("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("1.0"))
        ("r,repeat", "Repeat last run")
        ("h,help", "Print usage");

    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return 0;
    }

    unsigned normal_param;
    size_t events_param;
    float smooth_param;
    std::vector<Covariant<3>::Event> means;
    std::vector<float> stddevs;
    std::vector<float> fractions;

    if (result.count("repeat")) {
        std::cout << "Repeating previous run from file." << std::endl;
        std::ifstream infile("run_params.bin", std::ios::binary);
        if (!infile) {
            std::cerr << "Error: repeat file not found (run_params.bin). Run without --repeat first." << std::endl;
            return 1;
        }
        infile.read(reinterpret_cast<char*>(&normal_param), sizeof(normal_param));
        infile.read(reinterpret_cast<char*>(&events_param), sizeof(events_param));
        infile.read(reinterpret_cast<char*>(&smooth_param), sizeof(smooth_param));
        
        means.resize(normal_param);
        stddevs.resize(normal_param);
        if (normal_param > 1) {
            fractions.resize(normal_param - 1);
        }

        infile.read(reinterpret_cast<char*>(means.data()), means.size() * sizeof(Covariant<3>::Event));
        infile.read(reinterpret_cast<char*>(stddevs.data()), stddevs.size() * sizeof(float));
        if (!fractions.empty()) {
            infile.read(reinterpret_cast<char*>(fractions.data()), fractions.size() * sizeof(float));
        }
        infile.close();
    } else {
        normal_param = result["normal"].as<unsigned>();
        events_param = result["events"].as<size_t>();
        smooth_param = result["smooth"].as<float>();

        while (means.size() < normal_param) {
            Covariant<3>::Event mean {random_float() * .8f + .1f, random_float() * .8f + .1f, random_float() * .8f + .1f};
            means.push_back(mean);
            stddevs.push_back(random_float() * 0.05f + 0.025f);
        }
        if (normal_param > 1) {
            while (fractions.size() < normal_param - 1) {
                fractions.push_back(random_float());
            }
            std::sort(fractions.begin(), fractions.end());
        }

        std::ofstream outfile("run_params.bin", std::ios::binary | std::ios::trunc);
        outfile.write(reinterpret_cast<const char*>(&normal_param), sizeof(normal_param));
        outfile.write(reinterpret_cast<const char*>(&events_param), sizeof(events_param));
        outfile.write(reinterpret_cast<const char*>(&smooth_param), sizeof(smooth_param));
        outfile.write(reinterpret_cast<const char*>(means.data()), means.size() * sizeof(Covariant<3>::Event));
        outfile.write(reinterpret_cast<const char*>(stddevs.data()), stddevs.size() * sizeof(float));
        if (!fractions.empty()) {
            outfile.write(reinterpret_cast<const char*>(fractions.data()), fractions.size() * sizeof(float));
        }
        outfile.close();
    }
    
    std::cout << "Program running with --normal=" << normal_param << " and --events=" << events_param << " and --smooth=" << smooth_param << std::endl;
    
    Covariant<3> covariant(64);
    Covariant<3>::Event e;

    while (covariant.events() < events_param) {
        int population = std::upper_bound(fractions.begin(), fractions.end(), random_float()) - fractions.begin();
            for (unsigned i = 0; i < covariant.dimension; i++)
                e[i] = random_normal(means[population][i], stddevs[population]);
        covariant.event(e);
    }

    std::cout << "Analyzing the sample..." << std::endl;

    covariant.parameters(smooth_param);
    assert(covariant.factorProbability() < 0.01);
    assert(covariant.differentialEquation() < 0.01);

    std::cout << "Consistency checkes passed..." << std::endl;
    
    // Write the P array to a binary file
    std::cout << "Writing data to files..." << std::endl;
    auto write_file = [&](std::string filename, std::function<float(size_t)> get_val) {
        std::ofstream outfile(filename, std::ios::binary | std::ios::trunc);
        int dims[3] = {65, 65, 65};
        outfile.write(reinterpret_cast<char*>(&dims), sizeof(dims));
        for (size_t x = 0; x < covariant.size(); ++x) {
            float val = get_val(x);
            outfile.write(reinterpret_cast<char*>(&val), sizeof(float));
        }
        outfile.close();
    };

    write_file("P.bin", [&](size_t x) { return covariant.P(x); });
    write_file("M1.bin", [&](size_t x) { return covariant.M(0, x); });
    write_file("M2.bin", [&](size_t x) { return covariant.M(1, x); });
    write_file("M3.bin", [&](size_t x) { return covariant.M(2, x); });
    write_file("L.bin", [&](size_t x) { return covariant.L(x); });
    std::cout << "Data written successfully." << std::endl;
    
    return 0;
}