#include "Leonard.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

Marginal_Results Leonard::do_Marginal(std::shared_ptr<Laplace_Results> laplace, const std::vector<std::vector<float> *> &data, const Measurement i, const Measurement j, std::string pop_name) noexcept
{
    Marginal_Results results(laplace);

    Weighty<2> marginal(256);
    marginal.antialias = params.antialias;
    marginal.verify = params.verify;

    Coordinates<2> marginal_coord(marginal);
    Event<2> marginal_event;
    std::vector<unsigned short> marginal_klass(marginal.points(0) * marginal.points(1), 0);
    std::fill(marginal_klass.begin(), marginal_klass.end(), 0);

    for (unsigned short c = 0; c <= laplace->valid_clusters; ++c)
    {
        for (unsigned e : laplace->cluster_events[c])
        {
            marginal_event[0] = (*data[i])[laplace->idx[e]];
            marginal_event[1] = (*data[j])[laplace->idx[e]];
            marginal.event(marginal_event);
            if (c == 0)
                continue;

            marginal.locate(marginal_event, marginal_coord);
            size_t idx = (size_t)marginal_coord;
            unsigned short d = marginal_klass[idx];
            if (d == 0)
                marginal_klass[idx] = c;
            else
                marginal_klass[idx] = std::min(c, d);
        }
    }
    marginal.prepare(params.smoothing);

    auto class_data = std::make_shared<std::vector<std::vector<std::vector<double>>>>(3, std::vector<std::vector<double>>(marginal.points(0), std::vector<double>(marginal.points(1))));
    auto quant_data = std::make_shared<std::vector<std::vector<double>>>(marginal.points(0), std::vector<double>(marginal.points(1)));

    for (unsigned y = 0; y < marginal.points(1); ++y)
    {
        for (unsigned x = 0; x < marginal.points(0); ++x)
        {
            size_t idx = x + y * marginal.points(0);
            for (unsigned i = 0; i < 3; ++i)
                if (marginal_klass[idx] == 0)
                    (*class_data)[i][y][x] = 255;
                else
                {
                    unsigned hue = (255 * (marginal_klass[idx] - 1) / (laplace->valid_clusters - 1));
                    (*class_data)[i][y][x] = 255 * colors[hue][i];
                }
            (*quant_data)[y][x] = (double)static_cast<const Function<2, float> &>(marginal.quantile)[idx];
        }
    }

    std::string over_path = params.img_dir + "/sample_" + selections.variables[i] + "_" + selections.variables[j] + ".png";
    std::string under_path = params.img_dir + "/sample_" + selections.variables[i] + "_" + selections.variables[j] + "_under.png";
    write_rgb_png(under_path, *class_data);
    laplace->future_plots.push_back(plot_plane.enqueue([this, class_data, quant_data, over_path]()
                                                       { make_marginal_plot(over_path, *class_data, *quant_data); }));
    laplace->sample_images.push_back("images/sample_" + selections.variables[i] + "_" + selections.variables[j] + ".png");

    for (unsigned c = 1; c <= laplace->valid_clusters; ++c)
    {
        if (laplace->cluster_events[c].empty())
            continue;

        std::vector<unsigned short> marginal_klass(marginal.size(), 0);

        marginal.reset();
        for (auto &e : laplace->cluster_events[c])
        {
            marginal_event[0] = (*data[i])[laplace->idx[e]];
            marginal_event[1] = (*data[j])[laplace->idx[e]];
            marginal.event(marginal_event);
            if (c == 0)
                continue;

            marginal.locate(marginal_event, marginal_coord);
            size_t idx = (size_t)marginal_coord;
            marginal_klass[idx] = c;
        }
        marginal.prepare(selections.smoothing);

        // make cluster marginal plot
        auto class_data = std::make_shared<std::vector<std::vector<std::vector<double>>>>(3, std::vector<std::vector<double>>(marginal.points(0), std::vector<double>(marginal.points(1))));
        auto quant_data = std::make_shared<std::vector<std::vector<double>>>(marginal.points(0), std::vector<double>(marginal.points(1)));

        for (unsigned y = 0; y < marginal.points(1); ++y)
        {
            for (unsigned x = 0; x < marginal.points(0); ++x)
            {
                size_t idx = x + y * marginal.points(0);
                for (unsigned i = 0; i < 3; ++i)
                    if (marginal_klass[idx] == 0)
                        (*class_data)[i][y][x] = 255;
                    else
                    {
                        unsigned hue = (255 * (marginal_klass[idx] - 1) / (laplace->valid_clusters - 1));
                        (*class_data)[i][y][x] = 255 * colors[hue][i];
                    }
                (*quant_data)[y][x] = (double)static_cast<const Function<2, float> &>(marginal.quantile)[idx];
            }
        }

        std::string path = params.img_dir + "/cluster_" + std::to_string(c) + "_" + selections.variables[i] + "_" + selections.variables[j] + ".png";
        std::string under_path = params.img_dir + "/cluster_" + std::to_string(c) + "_" + selections.variables[i] + "_" + selections.variables[j] + "_under.png";
        write_rgb_png(under_path, *class_data);
        laplace->future_plots.push_back(plot_plane.enqueue([path, class_data, quant_data]()
                                                           { make_marginal_plot(path, *class_data, *quant_data); }));
        laplace->cluster_images.push_back("images/cluster_" + std::to_string(c) + "_" + selections.variables[i] + "_" + selections.variables[j] + ".png");
    }

    return results;
}
