#include <functional>
#include <fstream>
#include <filesystem>

#include <Eigen/Dense>
#include <cxxopts.hpp>

#include "Covariant.hpp"
#include "Samples.hpp"
#include "Gating.hpp"
#include "FlowJo.hpp"

struct Params
{
    std::vector<std::string> files;
    std::vector<std::string> variables;
    std::vector<std::string> populations;
    float smoothing = 0.01f;
    float threshold = 0.001f;
    unsigned max_clusters = 50;
    size_t min_events = 100;
    unsigned grid_size = 256;
    int analysis_choice = 1;
};

class Leonard
{
public:
    Params params;
    SelectionState selections;
    unsigned clusters_found = 0;
    size_t laplacian_offset = 0;
    Workspace ws;
    std::vector<SampleData> dummy_samples;

    int parse_args(int argc, char *argv[])
    {
        cxxopts::Options options("Leonard", "Laplacian and Riemannian analysis from FlowJo workspaces");

        options.add_options()
        ("f,files", "List of file names", cxxopts::value<std::vector<std::string>>())
        ("v,variables", "List of variables", cxxopts::value<std::vector<std::string>>())
        ("p,populations", "List of populations", cxxopts::value<std::vector<std::string>>())
        ("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("0.01"))
        ("t,threshold", "Threshold", cxxopts::value<float>()->default_value("0.001"))
        ("max-clusters", "Max clusters", cxxopts::value<unsigned>()->default_value("50"))
        ("min-events", "Min events", cxxopts::value<size_t>()->default_value("100"))
        ("g,grid", "Grid resolution", cxxopts::value<unsigned>()->default_value("256"))
        ("a,analysis", "Analysis choice (0=EPP, 1=Laplace)", cxxopts::value<int>()->default_value("0"))
        ("h,help", "Print usage");

        options.parse_positional({"files", "variables", "populations"});

        auto result = options.parse(argc, argv);

        if (result.count("help"))
        {
            std::cout << options.help() << std::endl;
            return 0;
        }

        if (result.count("files")) params.files = result["files"].as<std::vector<std::string>>();
        if (result.count("variables")) params.variables = result["variables"].as<std::vector<std::string>>();
        if (result.count("populations")) params.populations = result["populations"].as<std::vector<std::string>>();
        params.smoothing = result["smooth"].as<float>();
        params.threshold = result["threshold"].as<float>();
        params.max_clusters = result["max-clusters"].as<unsigned>();
        params.min_events = result["min-events"].as<size_t>();
        params.grid_size = result["grid"].as<unsigned>();
        params.analysis_choice = result["analysis"].as<int>();

        return -1;
    }

    template <unsigned Dimension>
    void do_it(const std::vector<std::vector<float> *> &data, const std::vector<bool> &included, std::vector<unsigned short> &classes)
    {
        std::vector<size_t> idx;
        idx.reserve(included.size());
        for (size_t i = 0; i < included.size(); ++i)
            if (included[i])
                idx.push_back(i);

        if (idx.empty()) {
            std::cout << "No valid events found." << std::endl;
            return;
        }

        Events<Dimension> events;
        events.resize(idx.size());
        for (unsigned d = 0; d < Dimension; ++d)
            for (size_t i = 0; i < events.size(); ++i)
                events[i][d] = (*data[d])[idx[i]];

        std::cout << "Processing " << events.size() << " events..." << std::endl;
        Laplace<Dimension> laplace(selections.grid_size);
        for (const auto &e : events)
            laplace.event(e);
        std::cout << "Found " << events.size() << " valid events..." << std::endl;

        std::cout << "Analyzing the sample..." << std::endl;
        laplace.analyze(selections.smoothing, selections.threshold);
        if (laplace.differentialError() > .0001)
            std::cout << "Differential equation solution is unusually bad " << laplace.differentialError() << std::endl;
        else
            std::cout << "Consistency checks passed..." << std::endl;

        std::cout << "Performing Laplacian clustering..." << std::endl;
        clusters_found = laplace.cluster(selections.threshold);
        std::fill(classes.begin(), classes.end(), 0);
        for (size_t i = 0; i < events.size(); ++i)
            classes[idx[i]] = laplace.classify(events[i]);
        std::cout << "Found " << clusters_found << " clusters." << std::endl;

        for (size_t i = 0; i < idx.size(); ++i)
            classes[idx[i]] = laplace.classify(events[i]) + laplacian_offset;
    }

    int run()
    {
        bool is_datafile = false;
        if (!params.files.empty() && params.files[0].find(".wsp") == std::string::npos)
            is_datafile = true;

        if (is_datafile)
        {
            if (params.variables.empty())
            {
                std::cerr << "Error: Variables parameter is required for datafiles.\n";
                return 1;
            }
            if (!params.populations.empty())
            {
                std::cerr << "Error: Populations is always All for datafiles.\n";
                return 1;
            }
            selections.populations.push_back("All");
            selections.variables = params.variables;
            selections.smoothing = params.smoothing;
            selections.threshold = params.threshold;
            selections.max_clusters = params.max_clusters;
            selections.min_events = params.min_events;
            selections.grid_size = params.grid_size;
            selections.analysis_choice = params.analysis_choice;

            for (const auto& file : params.files)
            {
                SampleData sd;
                sd.name = file;
                sd.selected = true;
                dummy_samples.push_back(sd);
            }
            for (auto& sd : dummy_samples)
            {
                selections.samples.push_back(&sd);
            }
        }
        else
        {
            std::string filename;
            if (params.files.empty())
            {
                filename = find_workspace(0, nullptr);
            }
            else
            {
                filename = params.files[0];
                if (!std::filesystem::exists(filename))
                {
                    if (std::filesystem::exists(filename + ".wsp"))
                        filename += ".wsp";
                }
            }

            if (filename.empty())
                return 1;

            ws = parse_workspace(filename);
            if (ws.samples.empty())
            {
                std::cerr << "Error: No samples found or failed to parse workspace.\n";
                return 1;
            }

            if (params.variables.empty() && params.populations.empty())
            {
                selections = build_ftxui_interface(ws);
                if (selections.cancelled)
                    return 1;
            }
            else
            {
                selections.variables = params.variables;
                selections.populations = params.populations;
                selections.smoothing = params.smoothing;
                selections.threshold = params.threshold;
                selections.max_clusters = params.max_clusters;
                selections.min_events = params.min_events;
                selections.grid_size = params.grid_size;
                selections.analysis_choice = params.analysis_choice;

                // the selected samples should be all the samples that have all the populations specified
                for (auto& s : ws.samples)
                {
                    bool has_all = true;
                    for (const auto& p : selections.populations)
                    {
                        if (std::find(s.populations.begin(), s.populations.end(), p) == s.populations.end())
                        {
                            has_all = false;
                            break;
                        }
                    }
                    if (has_all)
                    {
                        s.selected = true;
                        selections.samples.push_back(&s);
                    }
                }
            }
        }

        std::cout << "\n=== LEONARD SELECTIONS ===\n";
        if (!selections.samples.empty())
        {
            std::cout << "Samples:";
            for (const auto *s : selections.samples)
                std::cout << " " << s->name;
            std::cout << "\n";
        }
        std::cout << "Variables: ";
        for (const auto &v : selections.variables)
            std::cout << v << " ";

        std::cout << "\nPopulations: ";
        for (const auto &p : selections.populations)
            std::cout << p << " ";

        std::cout << "\n\nSettings:\n";
        std::cout << "  Smoothing: " << selections.smoothing << "\n";
        std::cout << "  Threshold: " << selections.threshold << "\n";
        std::cout << "  Max Clusters: " << selections.max_clusters << "\n";
        std::cout << "  Min Events: " << selections.min_events << "\n";
        std::cout << "  Grid Size: " << selections.grid_size << "\n";

        std::cout << "\nAnalysis Method:\n"
                  << analysis_choices[selections.analysis_choice] << "\n\n";

        for (const auto *s_ptr : selections.samples)
        {
            const auto &s = *s_ptr;
            std::unordered_map<std::string, std::shared_ptr<Gate>> gate_id_to_gate;
            for (const auto &g : s.gates)
                if (g)
                    gate_id_to_gate[g->id] = g;

            std::cout << "Processing sample: " << s.name << " ... ";
            DataSet dataset;
            if (!dataset.read(s.name)) {
                std::cerr << std::endl << "Failed to read dataset: " << s.name << std::endl;
                continue;
            }
            std::cout << dataset.size() << " events read" << std::endl;

            Classification &laplace = dataset.classification["Laplace"];
            if (!laplace.classifications)
                laplace.classifications = std::make_shared<std::vector<unsigned short>>(dataset.size());
            for (size_t i = 0; i < dataset.size(); ++i)
            { 
                (*laplace.classifications)[i] /= 16;
                if ((*laplace.classifications)[i] > laplacian_offset)
                    laplacian_offset = (*laplace.classifications)[i];
            }
            laplacian_offset++;

            // setup spectral unmixing
            std::vector<std::vector<float> *> detector_data;
            for (const auto &v : s.spillover_matrix.parameters)
                detector_data.push_back(dataset.variable[v].data.get());
            for (const auto &pair : s.spillover_matrix.comp_infix_map)
                Variable &var = dataset.variable[pair.second];

            if (!s.spillover_matrix.matrix.empty() && !s.spillover_matrix.matrix[0].empty()) {
                const std::vector<std::vector<double>> &spill = s.spillover_matrix.matrix;
                Eigen::MatrixXf spectrum(spill.size(), spill[0].size());
                for (size_t i = 0; i < spill.size(); i++)
                    for (size_t j = 0; j < spill[i].size(); j++)
                        spectrum(i, j) = spill[i][j];
                Eigen::BDCSVD<Eigen::MatrixXf, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(spectrum);
                Eigen::MatrixXf pseudo(svd.matrixV() * svd.singularValues().asDiagonal().inverse() * svd.matrixU().transpose());
                for (size_t i = 0; i < pseudo.rows(); i++)
                {
                    const auto &p = s.spillover_matrix.parameters[i];
                    auto v = s.spillover_matrix.comp_infix_map.find(p);
                    if (v == s.spillover_matrix.comp_infix_map.end())
                        continue;
                    Variable &var = dataset.variable[v->second];
                    for (size_t j = 0; j < pseudo.cols(); j++)
                        var.unmixing.push_back(pseudo(i, j));
                }
            }

            // set the scales on the variables for this sample
            for (auto &var : dataset.variable)
            {
                auto xfrm = s.transforms.find(var.first);
                if (xfrm != s.transforms.end())
                    var.second.transform = xfrm->second;
            }

            // cache of transformed data
            std::unordered_map<std::string, std::shared_ptr<std::vector<float>>> transformed_data;

            auto get_transformed = [&](const std::string &var_name) -> std::shared_ptr<std::vector<float>>
            {
                auto it = transformed_data.find(var_name);
                if (it != transformed_data.end())
                    return it->second;

                auto &var = dataset.variable[var_name];

                // if the data doesn't exist we need to do the unmixing first
                if (!var.data && !var.unmixing.empty())
                {
                    var.data = std::make_shared<std::vector<float>>(dataset.size(), 0.0f);
                    float *bare_data = var.data.get()->data();
                    for (size_t j = 0; j < detector_data.size(); ++j)
                    {
                        float *bare_detector = detector_data[j]->data();
                        float coefficient = var.unmixing[j];
                        for (size_t i = 0; i < dataset.size(); ++i)
                                bare_data[i] += bare_detector[i] * coefficient;
                    }
                }

                std::shared_ptr<std::vector<float>> transformed;
                const std::shared_ptr<std::vector<float>> &orig_data = var.data;
                if (!orig_data) {
                    return nullptr;
                }
                
                transformed = std::make_shared<std::vector<float>>(orig_data->size());
                if (var.transform) {
                    var.transform->transform(*orig_data, *transformed);
                } else {
                    *transformed = *orig_data;
                }

                transformed_data[var_name] = transformed; // Cache it for future lookups
                return transformed;
            };

            std::cout << "Scaling and unmixing data for analysis..." << std::endl;
            // generate the appropriate scaled data for analysis
            bool vars_valid = true;
            for (size_t v = 0; v < selections.variables.size(); ++v)
            {
                const auto &vname = selections.variables[v];
                Variable &var = dataset.variable[vname];
                auto transformed = get_transformed(vname);
                if (transformed)
                    var.data = transformed;
                else {
                    std::cerr << "Variable " << vname << " not found in sample " << s.name << std::endl;
                    vars_valid = false;
                }
            }
            if (!vars_valid) continue;

            std::cout << "Finding all the gates and, scaling and unmixing the needed data..." << std::endl;
            // find all the gates and any data they need
            std::vector<std::shared_ptr<Gate>> gates;
            for (size_t i = 0; i < selections.populations.size(); ++i)
            {
                std::string pop_name = selections.populations[i];
                if (pop_name == "All") {
                    Subset &subset = dataset.subset[pop_name];
                    subset.membership = std::make_shared<std::vector<bool>>(dataset.size(), true);
                    continue;
                }
                
                Subset &subset = dataset.subset[pop_name];
                std::string gate_id;
                std::shared_ptr<Gate> gate;
                auto it = s.pop_to_gate_id.find(pop_name);
                if (it != s.pop_to_gate_id.end())
                {
                    gate_id = it->second;
                    do
                    {
                        gate = gate_id_to_gate[gate_id];
                        auto it = std::find(gates.begin(), gates.end(), gate);
                        if (it != gates.end())
                            break;
                        gates.push_back(gate);
                        for (auto &dim : gate->dimensions)
                        {
                            dim.transform = dataset.variable[dim.name].transform;
                            gate->data.push_back(get_transformed(dim.name));
                        }
                        gate_id = gate->parent_id;
                    } while (!gate_id.empty());
                }
            }

            std::cout << "Evaluating all the necessary gates..." << std::endl;
            // evaluate all the gates needed
            for (auto it = gates.rbegin(); it != gates.rend(); ++it)
            {
                auto &gate = *it;
                Subset &subset = dataset.subset[gate->name];
                if (gate->parent_id.empty())
                    subset.membership = std::make_shared<std::vector<bool>>(dataset.size(), true);
                else
                {
                    auto parent_it = std::find(gates.begin(), gates.end(), gate_id_to_gate[gate->parent_id]);
                    Subset &parent = dataset.subset[gate_id_to_gate[gate->parent_id]->name];
                    subset.membership = std::make_shared<std::vector<bool>>(*parent.membership);
                }
                gate->evaluate(*subset.membership);
            }

            std::vector<unsigned short> *classification = nullptr;
            if (selections.analysis_choice == 1)
                classification = laplace.classifications.get();

            // scaled data marshalled as events for analysis
            int num_vars_selected = selections.variables.size();
            std::vector<std::vector<float> *> data(num_vars_selected);
            for (size_t i = 0; i < num_vars_selected; ++i)
                data[i] = dataset.variable[selections.variables[i]].data.get();

            for (std::string &pop_name : selections.populations)
            {
                std::cout << "Analysing population " << pop_name << std::endl;
                Subset &subset = dataset.subset[pop_name];
                std::vector<bool> *subpopulation = dataset.subset[pop_name].membership.get();
                if (!subpopulation) continue;

                switch (selections.analysis_choice)
                {
                case 0:
                    // EPP
                    break;
                case 1:
                    // Laplace
                    switch (num_vars_selected)
                    {
                    case 2:
                        do_it<2>(data, *subpopulation, *classification);
                        break;
                    case 3:
                        do_it<3>(data, *subpopulation, *classification);
                        break;
                    case 4:
                        do_it<4>(data, *subpopulation, *classification);
                        break;
                    default:
                        // error
                        break;
                    }
                    if (selections.analysis_choice == 1)
                    {
                        if (!is_datafile && !ws.filename.empty())
                            add_laplace_gates(ws.filename, s.name, pop_name, clusters_found, laplacian_offset);
                        laplacian_offset += clusters_found + 1;
                    }

                    break;
                default:
                    // error
                    break;
                }
            }

            if (selections.analysis_choice == 1)
            {
                std::filesystem::path len_path = s.name;
                std::string len_filename = len_path.stem().string();
                std::replace(len_filename.begin(), len_filename.end(), ' ', '_');
                len_path.replace_filename(len_filename + ".len");
                std::ofstream out(len_path);
                if (out.is_open())
                {
                    out << "Laplace\n";
                    for (auto val : *laplace.classifications)
                    {
                        out << 8 + 16 * val << "\n";
                    }
                    out.close();
                }

                if (!is_datafile && !ws.filename.empty())
                {
                    add_laplace_derived_parameter(ws.filename, s.name);
                }
            }
        };

        return 0;
    }
};

int main(int argc, char *argv[])
{
    Leonard leonard;
    int res = leonard.parse_args(argc, argv);
    if (res != -1) return res;

    return leonard.run();
}
