#include <iostream>
#include <string>
#include <fstream>
#include <functional>
#include <cxxopts.hpp>
#include "TestData.hpp"
#include "Samples.hpp"

struct Params
{
    std::string filename;
    std::vector<std::string> variables;
    std::vector<std::string> populations;
    unsigned dimension;
    unsigned normal;
    unsigned snake;
    unsigned exponential;
    size_t events;
    bool convert;
    bool ascii;
};

class TestDataApp
{
public:
    Params params;

    int parse_args(int argc, char *argv[])
    {
        // Set up the command-line options parser.
        cxxopts::Options options("WeightyCLI", "Generate test data for the Weighty class");

        // Add the command-line options.
        options.add_options()
        ("f,file", "File name for generated test data", cxxopts::value<std::string>()->default_value("test_data"))
        ("variables", "List of variables", cxxopts::value<std::vector<std::string>>())
        ("populations", "List of populations", cxxopts::value<std::vector<std::string>>())
        ("d,dimension", "Dimension of the events", cxxopts::value<unsigned>()->default_value("2"))
        ("n,normal", "Normal distributions", cxxopts::value<unsigned>()->default_value("3"))
        ("exponential", "Exponential distributions", cxxopts::value<unsigned>()->default_value("0"))
        ("snake", "Snake distributions", cxxopts::value<unsigned>()->default_value("0"))
        ("e,events", "Number of events to generate", cxxopts::value<size_t>()->default_value("10000"))
        ("convert", "Convert existing file between ASCII and binary", cxxopts::value<bool>()->default_value("false"))
        ("a,ascii", "Use ASCII format for data files", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
        ("h,help", "Print usage");

        options.parse_positional({"file", "variables", "populations"});

        auto result = options.parse(argc, argv);

        if (result.count("help"))
        {
            std::cout << options.help() << std::endl;
            return 0;
        }

        params.filename = result["file"].as<std::string>();
        if (result.count("variables")) {
            params.variables = result["variables"].as<std::vector<std::string>>();
            params.dimension = params.variables.size();
        } else {
            params.dimension = result["dimension"].as<unsigned>();
        }
        if (result.count("populations")) {
            params.populations = result["populations"].as<std::vector<std::string>>();
        }
        params.normal = result["normal"].as<unsigned>();
        params.snake = result["snake"].as<unsigned>();
        params.exponential = result["exponential"].as<unsigned>();
        params.events = result["events"].as<size_t>();
        params.convert = result["convert"].as<bool>();
        params.ascii = result["ascii"].as<bool>();

        return -1;
    }

    template <unsigned Dimension>
    int do_it()
{
    TestData<Dimension> test_data;
    if (params.convert)
    {
        bool source_is_ascii = !params.ascii;
        std::string source_ext = source_is_ascii ? ".txt" : ".dat";
        std::cout << "Converting " << params.filename << " to "
                  << (params.ascii ? "ASCII" : "binary") << "..." << std::endl;
        DataSet dataset;
        if (!dataset.read(params.filename))
        {
            std::cerr << "Error: Could not read source file for conversion: " << params.filename << std::endl;
            return 1;
        }
        std::vector<size_t> var_indices;
        if (!params.variables.empty()) {
            for (const auto& var_name : params.variables) {
                bool found = false;
                for (size_t i = 0; i < dataset.size(); ++i) {
                    if (dataset[i].name == var_name) {
                        var_indices.push_back(i);
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    std::cerr << "Error: Could not find variable '" << var_name << "' in dataset." << std::endl;
                    return 1;
                }
            }
        } else {
            for (size_t i = 0; i < params.dimension && i < dataset.size(); ++i) {
                var_indices.push_back(i);
            }
            if (var_indices.size() != params.dimension) {
                std::cerr << "Error: Dataset has fewer variables than requested dimension." << std::endl;
                return 1;
            }
        }

        std::vector<size_t> pop_indices;
        if (!params.populations.empty()) {
            for (const auto& pop_name : params.populations) {
                bool found = false;
                for (size_t i = 0; i < dataset.size(); ++i) {
                    if (dataset[i].name == pop_name) {
                        pop_indices.push_back(i);
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    std::cerr << "Error: Could not find population '" << pop_name << "' in dataset." << std::endl;
                }
            }
        }
        Projection proj(var_indices, dataset);
        size_t num_events = proj[0].get_data<float>().size();
        test_data.resize(num_events);
        for (unsigned d = 0; d < Dimension; ++d) {
            const auto& col = proj[d].get_data<float>();
            for (size_t i = 0; i < num_events; ++i) {
                test_data[i][d] = col[i];
            }
        }
        test_data.readTruth(params.filename + ".truth", source_is_ascii);
    }
    else
    {
        std::cout << "Program running with normal=" << params.normal << " snake=" << params.snake << " exponential=" << params.exponential << std::endl;
        std::cout << "Generating " << params.events << " events in " << params.dimension << " dimensions..." << std::endl;
        typename TestData<Dimension>::RandomSample test_sample;
        for (unsigned i = 0; i < params.normal; i++)
            test_sample.subpopulation(new typename TestData<Dimension>::Normal());
        for (unsigned i = 0; i < params.snake; i++)
            test_sample.subpopulation(new typename TestData<Dimension>::Snake());
        for (unsigned i = 0; i < params.exponential; i++)
            test_sample.subpopulation(new typename TestData<Dimension>::Exponential());
        test_data.generate(test_sample, params.events);
    }

    std::string target_ext = params.ascii ? ".txt" : ".dat";
    std::cout << "Saving " << test_data.size() << " events to " << params.filename << target_ext << "..." << std::endl;
    test_data.write(params.filename + target_ext, params.ascii);
    test_data.writeTruth(params.filename + ".truth", params.ascii);

    return 0;
}

    int run()
    {
    switch (params.dimension)
    {
    case 2:
        return do_it<2>();
    case 3:
        return do_it<3>();
    case 4:
        return do_it<4>();
    default:
        std::cerr << "Error: Unsupported dimension: " << params.dimension << std::endl;
        return 1;
    }
    }
};

int main(int argc, char *argv[])
{
    TestDataApp app;
    int res = app.parse_args(argc, argv);
    if (res != -1) return res;
    return app.run();
}