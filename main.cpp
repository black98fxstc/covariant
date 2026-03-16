#include <iostream>
#include <string>
#include <random>
#include <fstream>

// A popular, header-only library for command-line parsing.
// See setup instructions below.
#include <cxxopts.hpp>

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
    
    Covariant<3> covariant(64);
    Covariant<3>::Event e;
    std::vector<Covariant<3>::Event> means;
    std::vector<float> stddevs;
    while (means.size() < normal_param) {
        Covariant<3>::Event mean {random_float() * .8f + .1f, random_float() * .8f + .1f, random_float() * .8f + .1f};
        means.push_back(mean);
        stddevs.push_back(random_float() * 0.05f + 0.025f);
    }
    std::vector<float> fractions;
    while (fractions.size() < normal_param - 1) {
        fractions.push_back(random_float());
    }
    std::sort(fractions.begin(), fractions.end());
    
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

    std::cout << "Cinsistency checkes passed..." << std::endl;
    
    // Write the P array to a binary file
    std::cout << "Writing P data to file..." << std::endl;
    {
        std::ofstream outfile("P_data.bin", std::ios::binary | std::ios::trunc);
        int dims[3] = {65, 65, 65};
        outfile.write(reinterpret_cast<char*>(&dims), sizeof(dims));
        for(size_t x = 0; x < covariant.size(); ++x){
            float val = covariant.P(x);
            outfile.write(reinterpret_cast<char*>(&val), sizeof(float));
        }
        for(size_t x = 0; x < covariant.size(); ++x){
            float val = covariant.S(1,x);
            outfile.write(reinterpret_cast<char*>(&val), sizeof(float));
        }
        for(size_t x = 0; x < covariant.size(); ++x){
            float val = covariant.S(2,x);
            outfile.write(reinterpret_cast<char*>(&val), sizeof(float));
        }
        for(size_t x = 0; x < covariant.size(); ++x){
            float val = covariant.S(3,x);
            outfile.write(reinterpret_cast<char*>(&val), sizeof(float));
        }
        for(size_t x = 0; x < covariant.size(); ++x){
            float val = covariant.L(x);
            outfile.write(reinterpret_cast<char*>(&val), sizeof(float));
        }
        outfile.close();
    }
    std::cout << "Data written successfully." << std::endl;
    
    return 0;
}