#include <functional>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <memory>
#include <algorithm>
#include <future>
#include <atomic>
#include <thread>
#include <unordered_map>
#include <numbers>
#include <cmath>

#include <Eigen/Dense>
#include <fftw3.h>
#include <cxxopts.hpp>

#ifdef _WIN32
#include <windows.h>
#include <shellapi.h>
#endif

#include "Leonard.hpp"
#include "Reports.hpp"
#include "Samples.hpp"

int Leonard::parse_args(int argc, char *argv[])
{
    cxxopts::Options options("Leonard", "Laplacian and Riemannian analysis from FlowJo workspaces");

    options.add_options()("f,file", "File name", cxxopts::value<std::string>())("v,variables", "List of variables", cxxopts::value<std::string>())("p,populations", "List of populations", cxxopts::value<std::string>())("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("0.01"))("t,threshold", "Threshold", cxxopts::value<float>()->default_value("0.001"))("max-clusters", "Max clusters", cxxopts::value<unsigned>()->default_value("12"))("min-events", "Min cluster abs", cxxopts::value<size_t>()->default_value("0"))("min-relative", "Min cluster rel", cxxopts::value<float>()->default_value("0.0"))("kld-norm", "KLD Normal", cxxopts::value<float>()->default_value("0.04"))("kld-exp", "KLD Exponential", cxxopts::value<float>()->default_value("0.2"))("tolerance", "Tolerance", cxxopts::value<float>()->default_value("0.01"))("antialias", "Antialiasing", cxxopts::value<bool>()->default_value("true")->implicit_value("true"))("verify", "Verify consistency", cxxopts::value<bool>()->default_value("true")->implicit_value("true"))("g,grid", "Grid resolution", cxxopts::value<unsigned>()->default_value("256"))("a,analysis", "Analysis choice (0=EPP, 1=Laplace)", cxxopts::value<int>()->default_value("0"))("quiet", "Do not echo conversation to the console", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))("verbose", "Also echo the detailed log to the console", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))("h,help", "Print usage");

    options.parse_positional({"file", "variables", "populations"});

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        return 0;
    }

    if (result.count("file"))
        params.files.push_back(result["file"].as<std::string>());
    if (result.count("variables"))
    {
        std::stringstream ss(result["variables"].as<std::string>());
        std::string item;
        while (std::getline(ss, item, ','))
            params.variables.push_back(item);
    }
    if (result.count("populations"))
    {
        std::stringstream ss(result["populations"].as<std::string>());
        std::string item;
        while (std::getline(ss, item, ','))
            params.populations.push_back(item);
    }
    params.smoothing = result["smooth"].as<float>();
    params.threshold = result["threshold"].as<float>();
    params.max_clusters = result["max-clusters"].as<unsigned>();
    params.min_events = result["min-events"].as<size_t>();
    params.min_relative = result["min-relative"].as<float>();
    params.kld_normal = result["kld-norm"].as<float>();
    params.kld_exponential = result["kld-exp"].as<float>();
    params.tolerance = result["tolerance"].as<float>();
    params.grid_size = result["grid"].as<unsigned>();
    params.antialias = result["antialias"].as<bool>();
    params.verify = result["verify"].as<bool>();
    params.analysis_choice = result["analysis"].as<int>();
    session_.set_quiet(result["quiet"].as<bool>());
    session_.set_verbose(result["verbose"].as<bool>());

    return -1;
};

/**
 * @brief Simple HSV to RGB conversion helper.
 *
 * @param h Hue [0, 1]
 * @param s Saturation [0, 1]
 * @param v Value [0, 1]
 * @return std::vector<double> RGB components
 */

static std::vector<double> hsv_to_rgb(double h, double s, double v)
{
    double r = 0, g = 0, b = 0;
    if (s == 0)
    {
        r = g = b = v;
    }
    else
    {
        double h_pos = (h == 1.0) ? 0.0 : h * 6.0;
        int i = static_cast<int>(std::floor(h_pos));
        double f = h_pos - i;
        double p = v * (1.0 - s);
        double q = v * (1.0 - (s * f));
        double t = v * (1.0 - (s * (1.0 - f)));
        switch (i)
        {
        case 0:
            r = v;
            g = t;
            b = p;
            break;
        case 1:
            r = q;
            g = v;
            b = p;
            break;
        case 2:
            r = p;
            g = v;
            b = t;
            break;
        case 3:
            r = p;
            g = q;
            b = v;
            break;
        case 4:
            r = t;
            g = p;
            b = v;
            break;
        default:
            r = v;
            g = p;
            b = q;
            break;
        }
    }
    return {r, g, b};
}

int Leonard::run()
{
    configure_noninteractive_gnuplot();

    std::string wisdom_path;
#ifdef _WIN32
    if (const char *appdata = std::getenv("APPDATA"))
    {
        wisdom_path = std::string(appdata) + "\\Leonard";
    }
#else
    if (const char *home = std::getenv("HOME"))
    {
        wisdom_path = std::string(home) + "/.config/Leonard";
    }
#endif
    if (!wisdom_path.empty())
    {
        std::error_code ec;
        std::filesystem::create_directories(wisdom_path, ec);
        wisdom_path += "/wisdom.fftwf";
        if (std::filesystem::exists(wisdom_path))
        {
            fftwf_import_wisdom_from_filename(wisdom_path.c_str());
        }
    }

    if (fftwf_init_threads())
        fftwf_plan_with_nthreads(std::max(1u, std::thread::hardware_concurrency()));

    struct WisdomSaver
    {
        std::string path;
        ~WisdomSaver()
        {
            if (!path.empty())
            {
                std::lock_guard<std::mutex> lock(get_fftw_mutex());
                fftwf_export_wisdom_to_filename(path.c_str());
            }
        }
    } saver{wisdom_path};

    // create color map
    colors.reserve(256);
    for (size_t i = 0; i < 256; ++i)
    {
        double h = 0.6666 * (1.0 - static_cast<double>(i) / 256);
        colors.push_back(hsv_to_rgb(h, 1.0, 0.75));
    }

    bool is_datafile = false;
    if (!params.files.empty() && params.files[0].find(".wsp") == std::string::npos)
        is_datafile = true;

    if (is_datafile)
    {
        if (!params.files.empty())
        {
            std::filesystem::path df_path = std::filesystem::absolute(params.files[0]);
            if (df_path.has_parent_path())
            {
                std::error_code ec;
                std::filesystem::current_path(df_path.parent_path(), ec);
            }
        }

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
        selections.min_cluster_rel = params.min_relative;
        selections.kld_norm = params.kld_normal;
        selections.kld_exp = params.kld_exponential;
        selections.tolerance = params.tolerance;
        selections.grid_size = params.grid_size;
        selections.analysis_choice = params.analysis_choice;

        for (const auto &file : params.files)
        {
            SampleData sd;
            sd.name = file;
            sd.selected = true;
            dummy_samples.push_back(sd);
        }
        for (auto &sd : dummy_samples)
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

        std::filesystem::path wsp_path = std::filesystem::absolute(filename);
        if (wsp_path.has_parent_path())
        {
            std::error_code ec;
            std::filesystem::current_path(wsp_path.parent_path(), ec);
        }

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
            selections.min_cluster_rel = params.min_relative;
            selections.kld_norm = params.kld_normal;
            selections.kld_exp = params.kld_exponential;
            selections.tolerance = params.tolerance;
            selections.grid_size = params.grid_size;
            selections.analysis_choice = params.analysis_choice;

            // the selected samples should be all the samples that have all the populations specified
            for (auto &s : ws.samples)
            {
                bool has_all = true;
                for (const auto &p : selections.populations)
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

    std::filesystem::path source(!ws.filename.empty() ? ws.filename : (params.files.empty() ? "unknown" : params.files[0]));
    std::string report_dir = source.stem().string() + ".len";
    if (!session_.open(report_dir))
        std::cerr << "Could not open session log in " << report_dir << std::endl;
    else
        say << "Session log: " << session_.path().string() << std::endl;

    {
        auto message = say << "\n=== LEONARD SELECTIONS ===\n";
        if (!selections.samples.empty())
        {
            message << "Samples:";
            for (const auto *s : selections.samples)
                message << " " << s->name;
            message << "\n";
        }
        message << "Variables: ";
        for (const auto &v : selections.variables)
            message << v << " ";

        message << "\nPopulations: ";
        for (const auto &p : selections.populations)
            message << p << " ";

        message << "\n\nSettings:\n";
        message << "  Smoothing: " << selections.smoothing << "\n";
        message << "  Threshold: " << selections.threshold << "\n";
        message << "  Max Clusters: " << selections.max_clusters << "\n";
        message << "  Min Cluster Abs: " << selections.min_events << "\n";
        message << "  Min Cluster Rel: " << selections.min_cluster_rel << "\n";
        message << "  KLD Norm: " << selections.kld_norm << "\n";
        message << "  KLD Exp: " << selections.kld_exp << "\n";
        message << "  Tolerance: " << selections.tolerance << "\n";
        message << "  Grid Size: " << selections.grid_size << "\n";

        message << "\nAnalysis Method:\n"
                << analysis_choices[selections.analysis_choice] << "\n\n";
    }

    params.img_dir = std::filesystem::absolute(std::filesystem::path(report_dir) / "images")
                         .lexically_normal()
                         .generic_string();
    std::filesystem::create_directories(params.img_dir);
    std::vector<Reports::ReportLink> report_links;

    for (const auto *s_ptr : selections.samples)
    {
        const auto &s = *s_ptr;
        std::unordered_map<std::string, std::shared_ptr<Gate>> gate_id_to_gate;
        for (const auto &g : s.gates)
            if (g)
                gate_id_to_gate[g->id] = g;

        say << "Processing sample: " << s.name << " ... " << std::endl;
        DataSet dataset;
        if (!dataset.read(s.name))
        {
            session_.error("Failed to read dataset: " + s.name + "\n");
            continue;
        }
        say << dataset.size() << " events read" << std::endl;

        Classification &laplace = dataset.classification["Laplace"];
        if (!laplace.classifications)
            laplace.classifications = std::make_shared<std::vector<unsigned short>>(dataset.size());
        for (size_t i = 0; i < dataset.size(); ++i)
        {
            (*laplace.classifications)[i] /= 4;
            if ((*laplace.classifications)[i] > laplacian_offset)
                laplacian_offset = (*laplace.classifications)[i];
        }
        if (laplacian_offset > 0)
            laplacian_offset++;

        // setup spectral unmixing
        std::vector<std::vector<float> *> detector_data;
        for (const auto &v : s.spillover_matrix.parameters)
            detector_data.push_back(dataset.variable[v].data.get());
        for (const auto &pair : s.spillover_matrix.comp_infix_map)
            Variable &var = dataset.variable[pair.second];

        if (!s.spillover_matrix.matrix.empty() && !s.spillover_matrix.matrix[0].empty())
        {
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
                Eigen::Map<Eigen::ArrayXf> dest_map(var.data->data(), dataset.size());
                for (size_t j = 0; j < detector_data.size(); ++j)
                {
                    Eigen::Map<const Eigen::ArrayXf> src_map(detector_data[j]->data(), dataset.size());
                    dest_map += src_map * var.unmixing[j];
                }
            }

            std::shared_ptr<std::vector<float>> transformed;
            const std::shared_ptr<std::vector<float>> &orig_data = var.data;
            if (!orig_data)
            {
                return nullptr;
            }

            transformed = std::make_shared<std::vector<float>>(orig_data->size());
            if (var.transform)
            {
                var.transform->transform(*orig_data, *transformed);
            }
            else
            {
                *transformed = *orig_data;
            }

            transformed_data[var_name] = transformed; // Cache it for future lookups
            return transformed;
        };

        say << "Scaling and unmixing data for analysis..." << std::endl;
        // generate the appropriate scaled data for analysis
        bool vars_valid = true;
        for (size_t v = 0; v < selections.variables.size(); ++v)
        {
            const auto &vname = selections.variables[v];
            Variable &var = dataset.variable[vname];
            auto transformed = get_transformed(vname);
            if (transformed)
                var.data = transformed;
            else
            {
                session_.error("Variable " + vname + " not found in sample " + s.name + "\n");
                vars_valid = false;
            }
        }
        if (!vars_valid)
            continue;

        say << "Finding all the gates and, scaling and unmixing the needed data..." << std::endl;
        // find all the gates and any data they need
        std::vector<std::shared_ptr<Gate>> gates;
        for (size_t i = 0; i < selections.populations.size(); ++i)
        {
            std::string pop_name = selections.populations[i];
            if (pop_name == "All")
            {
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

        say << "Evaluating all the necessary gates..." << std::endl;
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

        std::vector<std::pair<std::string, std::future<Pursuit_Results>>> epp_results;
        std::vector<std::pair<std::string, std::future<std::shared_ptr<Laplace_Results>>>> laplace_results;

        for (std::string &pop_name : selections.populations)
        {
            say << "Enqueuing population " << pop_name << std::endl;
            std::vector<bool> subpopulation = *(dataset.subset[pop_name].membership.get());

            if (selections.analysis_choice == 0)
            {
                for (unsigned i = 0; i < num_vars_selected; i++)
                    for (auto it = std::find(subpopulation.begin(), subpopulation.end(), true); it != subpopulation.end(); it = std::find(it + 1, subpopulation.end(), true))
                        if ((*data[i])[it - subpopulation.begin()] < 0.0f || (*data[i])[it - subpopulation.begin()] > 1.0f)
                            subpopulation[it - subpopulation.begin()] = false;
                size_t total_events = std::count(subpopulation.begin(), subpopulation.end(), true);
                epp_results.push_back({pop_name, control_plane.enqueue([this, data, subpop = std::move(subpopulation), pop_name, total_events]() mutable
                                                                       { return do_Pursuit(data, std::move(subpop), pop_name, total_events); })});
            }
            else if (selections.analysis_choice == 1)
            {
                switch (num_vars_selected)
                {
                case 2:
                    laplace_results.push_back({pop_name, compute_plane.enqueue([this, data, subpopulation, pop_name]()
                                                                               { return do_Laplace<2>(data, subpopulation, pop_name); })});
                    break;
                case 3:
                    laplace_results.push_back({pop_name, compute_plane.enqueue([this, data, subpopulation, pop_name]()
                                                                               { return do_Laplace<3>(data, subpopulation, pop_name); })});
                    break;
                case 4:
                    laplace_results.push_back({pop_name, compute_plane.enqueue([this, data, subpopulation, pop_name]()
                                                                               { return do_Laplace<4>(data, subpopulation, pop_name); })});
                    break;
                }
            }
        }

        for (auto &result_pair : epp_results)
        {
            Pursuit_Results res = result_pair.second.get();
            res.wait_for_results();

            std::string filename = Reports::generate_epp_report(report_dir, s.name, result_pair.first, res, selections.variables);
            report_links.push_back({s.name + " - " + result_pair.first + " (EPP)", filename, "EPP tree analysis"});

            res.wait_for_plots();
        }
        for (auto &result_pair : laplace_results)
        {
            std::shared_ptr<Laplace_Results> res = result_pair.second.get();
            res->wait_for_results();
            if (res->idx.empty())
                continue;

            // Merge local classes into the global classification tracking synchronously
            for (size_t i = 0; i < res->idx.size(); ++i)
            {
                (*laplace.classifications)[res->idx[i]] = res->classification[i] + laplacian_offset;
            }

            if (!is_datafile && !ws.filename.empty())
            {
                add_laplace_gates(ws.filename, s.name, res->parent_name, res->clusters_found, laplacian_offset, res->cluster_events, selections.variables);
            }

            laplacian_offset += res->clusters_found + 1;
            if (res->cluster_events[res->valid_clusters + 1].size() > 0)
                ++laplacian_offset;

            res->wait_for_plots();

            for (const auto &entry : std::filesystem::directory_iterator(params.img_dir))
            {
                if (!entry.is_regular_file() || entry.path().extension() != ".png")
                    continue;
                if (entry.path().filename().string().find("_under.png") != std::string::npos ||
                    entry.path().filename().string().find(".rgba.png") != std::string::npos)
                    continue;
                make_overlay_transparent(entry.path().string());
            }

            std::string filename = Reports::generate_laplace_report(report_dir, s.name, result_pair.first, *res, selections.variables);
            report_links.push_back({s.name + " - " + result_pair.first + " (Laplace)", filename, "Laplacian clustering analysis"});
        }

        if (selections.analysis_choice == 1)
        {
            std::filesystem::path csv_path = s.name;
            std::string csv_filename = csv_path.stem().string();
            std::replace(csv_filename.begin(), csv_filename.end(), ' ', '_');
            csv_path.replace_filename(csv_filename + ".csv");
            std::ofstream out(csv_path);
            if (out.is_open())
            {
                out << "Laplace\n";
                for (auto val : *laplace.classifications)
                {
                    out << 2 + 4 * val << "\n";
                }
                out.close();
            }

            if (!is_datafile && !ws.filename.empty())
            {
                add_laplace_derived_parameter(ws.filename, s.name);
            }
        }
    };

    if (!report_links.empty())
    {
        Reports::update_index(report_dir, "Leonard Analysis Report", report_links);
        std::filesystem::path index_path = std::filesystem::absolute(std::filesystem::path(report_dir) / "index.html");
#if defined(_WIN32)
        ShellExecuteA(NULL, "open", index_path.string().c_str(), NULL, NULL, SW_SHOWNORMAL);
#elif defined(__APPLE__)
        std::string cmd = "open \"" + index_path.string() + "\"";
        std::system(cmd.c_str());
#else
        std::string cmd = "xdg-open \"" + index_path.string() + "\"";
        std::system(cmd.c_str());
#endif
    }

    return 0;
}

int main(int argc, char *argv[])
{
    Leonard leonard;
    int res = leonard.parse_args(argc, argv);
    if (res != -1)
        return res;

    return leonard.run();
}
