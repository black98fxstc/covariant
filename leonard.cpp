#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cstdlib>
#include <filesystem>
#include <algorithm>
#include <map>
#include <cctype>
#include <functional>

#include <Eigen/Dense>

#include "Covariant.hpp"
#include "Samples.hpp"
#include "Gating.hpp"
#include "FlowJo.hpp"

/*
Make leonard a command line program along the lines of darwin
craate a Leonard class to hold state including commandline Params
Three positional arguments a file name, a list of variables, and a list of populations all optional
If nothing at all is given look for a workspace file
if a workspace file was given or found construct and launch the tui
if a workpace file is given with arguments parse them into selections
the selected samples should be all the samples that have all the populations specified
if they start with a datafile the the second parameter is required and the third forbidden
and should be All in the selections.
move the main code into a method of Leonard and call it from main when the selections are ready to go
*/


template <unsigned Dimension>
void do_it(const std::vector<std::vector<float> *> &data, const std::vector<bool> &included, std::vector<unsigned short> &classes, SelectionState &params)
{
    std::vector<size_t> idx;
    for (size_t i = 0; i < included.size(); ++i)
        if (included[i])
            idx.push_back(i);

    Events<Dimension> events;
    events.resize(idx.size());
    for (unsigned d = 0; d < Dimension; ++d)
        for (size_t i = 0; i < events.size(); ++i)
            events[i][d] = (*data[d])[idx[i]];

    std::cout << "Processing " << events.size() << " events..." << std::endl;
    Laplace<Dimension> laplace(params.grid_size);
    for (const auto &e : events)
        laplace.event(e);
    std::cout << "Found " << events.size() << " valid events..." << std::endl;

    std::cout << "Analyzing the sample..." << std::endl;
    laplace.analyze(params.smoothing, params.threshold);
    if (laplace.differentialError() > .0001)
        std::cout << "Differential equation solution is unusually bad " << laplace.differentialError() << std::endl;
    else
        std::cout << "Consistency checks passed..." << std::endl;

    std::cout << "Performing Laplacian clustering..." << std::endl;
    unsigned found = laplace.cluster(params.threshold);
    std::fill(classes.begin(), classes.end(), 0);
    for (size_t i = 0; i < events.size(); ++i)
        classes[idx[i]] = laplace.classify(events[i]);
    std::cout << "Found " << found << " clusters." << std::endl;
};


int main(int argc, char *argv[])
{
    std::string filename = find_workspace(argc, argv);
    if (filename.empty())
        return 1;

    Workspace ws = parse_workspace(filename);
    if (ws.samples.empty())
    {
        std::cerr << "Error: No samples found or failed to parse workspace.\n";
        return 1;
    }

    SelectionState selections = build_ftxui_interface(ws);
    if (selections.cancelled)
        return 1;
    DataSet dataset;
    for (const auto *s_ptr : selections.samples)
    {
        const auto &s = *s_ptr;
        std::unordered_map<std::string, std::shared_ptr<Gate>> gate_id_to_gate;
        for (const auto &g : s.gates)
            if (g)
                gate_id_to_gate[g->id] = g;

        dataset.read(s.name);

        // setup spectral unmixing
        std::vector<std::vector<float> *> detector_data;
        for (const auto &v : s.spillover_matrix.parameters)
            detector_data.push_back(dataset.variable[v].data.get());
        for (const auto &pair : s.spillover_matrix.comp_infix_map)
            Variable &var = dataset.variable[pair.second];

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

            // if the data don't exist we need to do the unmixing first
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
            transformed = std::make_shared<std::vector<float>>(orig_data->size());
            var.transform->transform(*orig_data, *transformed);

            transformed_data[var_name] = transformed; // Cache it for future lookups
            return transformed;
        };

        // generate the appropriate scaled data for analysis
        for (size_t v = 0; v < selections.variables.size(); ++v)
        {
            const auto &vname = selections.variables[v];
            Variable &var = dataset.variable[vname];
            auto transformed = get_transformed(vname);
            if (transformed)
                var.data = transformed;
        }

        // find all the gates and any data they need
        std::vector<std::shared_ptr<Gate>> gates;
        for (size_t i = 0; i < selections.populations.size(); ++i)
        {
            std::string pop_name = selections.populations[i];
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

        // evaluate all the gates needed
        std::vector<bool> whole(dataset.size(), true);
        for (auto it = gates.rbegin(); it != gates.rend(); ++it)
        {
            auto &gate = *it;
            Subset &subset = dataset.subset[gate->name];
            if (gate->parent_id.empty())
                subset.membership = std::make_shared<std::vector<bool>>(dataset.size(), true);
            else
            {
                auto it = std::find(gates.begin(), gates.end(), gate_id_to_gate[gate->parent_id]);
                Subset &parent = dataset.subset[gate_id_to_gate[gate->parent_id]->name];
                subset.membership = std::make_shared<std::vector<bool>>(*parent.membership);
            }
            gate->evaluate(*subset.membership);
        }

        std::vector<unsigned short> *classification = nullptr;
        if (selections.analysis_choice == 1)
        {
            // someplace for the data to go
            Classification &laplace = dataset.classification["Laplace"];
            laplace.classifications = std::make_shared<std::vector<unsigned short>>(dataset.size(), 0);
            classification = laplace.classifications.get();
        }

        // scaled data marshalled as events for analysis
        int num_vars_selected = selections.variables.size();
        std::vector<std::vector<float> *> data(num_vars_selected);
        for (size_t i = 0; i < num_vars_selected; ++i)
            data[i] = dataset.variable[selections.variables[i]].data.get();

        for (std::string &pop_name : selections.populations)
        {
            std::vector<bool> *subpopulation = dataset.subset[pop_name].membership.get();
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
                    do_it<2>(data, *subpopulation, *classification, selections);
                    break;
                case 3:
                    do_it<3>(data, *subpopulation, *classification, selections);
                    break;
                case 4:
                    do_it<4>(data, *subpopulation, *classification, selections);
                    break;
                default:
                    // error
                    break;
                }
                break;
            default:
                // error
                break;
            }
        }
    };

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

    std::cout << "\nAnalysis Method:\n"
              << analysis_choices[selections.analysis_choice] << "\n\n";
    return 0;
}
