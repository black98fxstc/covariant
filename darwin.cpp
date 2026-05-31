#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cxxopts.hpp>
#include <filesystem>
#include <iomanip>
#include <nlohmann/json.hpp>
#include <cstdlib>
#include <fcntl.h>
#include <matplot/matplot.h>
#include <thread>
#include <algorithm>
#include <sstream>
#include <cctype>
#include <tuple>
#include <cmath>
#ifndef _WIN32
    #include <unistd.h>
    #include <sys/wait.h>
#endif

#include <libxml/parser.h>
#include <libxslt/xslt.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

#include "Covariant.hpp"

using json = nlohmann::json;

namespace fs = std::filesystem;

struct Params
{
    std::string filename;
    std::string out_dir;
    std::string log_path;
    std::string img_dir;
    unsigned dimension;
    unsigned grid;
    unsigned points;
    float smooth;
    float threshold;
    bool visual;
    unsigned max_clusters;
    size_t min_events;
    bool verbose;
    bool antialias;
    bool verify;
    bool grow;
    bool ascii;
    unsigned max_concurrent = 1;
    std::string labels_str;
    std::string style;
};

template <unsigned Dimension>
void for_each_plane(std::function<void(const unsigned i, const unsigned j)> func)
{
    for (unsigned i = 0; i < Dimension - 1; i++)
        for (unsigned j = i + 1; j < Dimension; j++)
            func(i, j);
};

void apply_stylesheet(const std::string& xml_path, const std::string& xsl_path, const std::string& out_path) {
    xmlDocPtr xml_doc = xmlParseFile(xml_path.c_str());
    if (!xml_doc) {
        std::cerr << "Error parsing XML file: " << xml_path << std::endl;
        return;
    }
    xsltStylesheetPtr xsl_doc = xsltParseStylesheetFile((const xmlChar*)xsl_path.c_str());
    if (!xsl_doc) {
        std::cerr << "Error parsing XSLT file: " << xsl_path << std::endl;
        xmlFreeDoc(xml_doc);
        return;
    }
    xmlDocPtr res_doc = xsltApplyStylesheet(xsl_doc, xml_doc, NULL);
    if (res_doc) {
        if (xsltSaveResultToFilename(out_path.c_str(), res_doc, xsl_doc, 0) == -1) {
            std::cerr << "Error writing output file: " << out_path << std::endl;
        }
        xmlFreeDoc(res_doc);
    } else {
        std::cerr << "Error applying stylesheet." << std::endl;
    }
    xsltFreeStylesheet(xsl_doc);
    xmlFreeDoc(xml_doc);
}

class Darwin
{
public:
    Params params;
    unsigned num_clusters = 0;
    std::ofstream xml_out;
    json report;
    std::vector<std::string> labels;
    std::vector<std::vector<double>> colors;
#ifndef _WIN32
    std::vector<pid_t> active_pids;
#endif

    Darwin(const Params& p) : params(p) {}

    ~Darwin() {
        wait_all();
    }
    /**
     * @brief Simple HSV to RGB conversion helper.
     * 
     * @param h Hue [0, 1]
     * @param s Saturation [0, 1]
     * @param v Value [0, 1]
     * @return std::vector<double> RGB components
     */
    static std::vector<double> hsv_to_rgb(double h, double s, double v) {
        double r = 0, g = 0, b = 0;
        if (s == 0) {
            r = g = b = v;
        } else {
            double h_pos = (h == 1.0) ? 0.0 : h * 6.0;
            int i = static_cast<int>(std::floor(h_pos));
            double f = h_pos - i;
            double p = v * (1.0 - s);
            double q = v * (1.0 - (s * f));
            double t = v * (1.0 - (s * (1.0 - f)));
            switch (i) {
                case 0: r = v; g = t; b = p; break;
                case 1: r = q; g = v; b = p; break;
                case 2: r = p; g = v; b = t; break;
                case 3: r = p; g = q; b = v; break;
                case 4: r = t; g = p; b = v; break;
                default: r = v; g = p; b = q; break;
            }
        }
        return {r, g, b};
    }

    void setup() {
        std::string header_line = params.labels_str;
        if (header_line.empty() && params.ascii) {
            std::string ext = params.ascii ? ".txt" : ".dat";
            std::ifstream in(params.filename + ext);
            if (in) {
                while (in.peek() != EOF && std::isspace(static_cast<unsigned char>(in.peek()))) {
                    in.ignore();
                }
                int first = in.peek();
                if (first != EOF && !std::isdigit(static_cast<unsigned char>(first)) && first != '-' && first != '+' && first != '.') {
                    std::getline(in, header_line);
                }
            }
        }
        
        if (!header_line.empty()) {
            std::replace(header_line.begin(), header_line.end(), ',', ' ');
            std::istringstream iss(header_line);
            std::string label;
            while (iss >> label) {
                std::string sanitized = label;
                for (char &c : sanitized) {
                    if (!std::isalnum(static_cast<unsigned char>(c))) c = '_';
                }
                if (!sanitized.empty()) {
                    labels.push_back(sanitized);
                }
            }
        }

        std::vector<std::string> default_labels = {"X", "Y", "Z", "W"};
        for (size_t i = labels.size(); i < params.dimension; ++i) {
            if (i < default_labels.size()) {
                labels.push_back(default_labels[i]);
            } else {
                labels.push_back("Dim" + std::to_string(i));
            }
        }
        labels.resize(params.dimension);

        // create color map
        colors.reserve(256);
        for (size_t i = 0; i < 256; ++i) 
        {
            double h = 0.6666 * (1.0 - static_cast<double>(i) / 256);
            colors.push_back(hsv_to_rgb(h, 1.0, 0.75));
        }

        // Create output directories
        fs::create_directories(params.out_dir);
        fs::create_directories(params.img_dir);

        // Create logs folder and truncate the session log file
        std::string log_dir = params.out_dir + "/logs";
        fs::create_directories(log_dir);
        params.log_path = log_dir + "/gnuplot.log";
        std::ofstream(params.log_path, std::ios::trunc);

        report["filename"] = params.filename;
        report["dimension"] = params.dimension;
        report["grid_size"] = params.grid;
        report["labels"] = labels;
        report["clusters"] = json::array();

        xml_out.open(params.out_dir + "/report.xml");
        xml_out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        xml_out << "<?xml-stylesheet type=\"text/xsl\" href=\"report.xsl\"?>\n";
        xml_out << "<DarwinReport file=\"" << params.filename << "\" dimensions=\"" << params.dimension << "\">\n";

        // Generate the XSLT stylesheet right next to the XML file
        std::ofstream xslt_out(params.out_dir + "/report.xsl");
        xslt_out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                 << "<xsl:stylesheet version=\"1.0\" xmlns:xsl=\"http://www.w3.org/1999/XSL/Transform\">\n"
                 << "  <xsl:output method=\"html\" indent=\"yes\"/>\n"
                 << "  <xsl:template match=\"/DarwinReport\">\n"
                 << "    <html>\n"
                 << "      <head>\n"
                 << "        <title>Darwin Report: <xsl:value-of select=\"@file\"/></title>\n"
                 << "        <style>\n"
                 << "          body { font-family: sans-serif; margin: 20px; background: #f4f4f9; color: #333; }\n"
                 << "          .row { display: flex; flex-direction: row; background: #fff; margin-bottom: 20px; padding: 15px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }\n"
                 << "          .info { flex: 0 0 350px; padding-right: 20px; border-right: 1px solid #eee; margin-right: 20px; }\n"
                 << "          .plots { display: flex; flex-wrap: wrap; gap: 10px; align-items: center; }\n"
                 << "          .plot img { width: 300px; height: auto; border: 1px solid #ddd; border-radius: 4px; }\n"
                 << "          .data-table { border-collapse: collapse; margin-top: 10px; width: 100%; font-size: 0.9em; }\n"
                 << "          .data-table th, .data-table td { border: 1px solid #ccc; padding: 6px; text-align: center; }\n"
                 << "          .data-table th { background: #f9f9f9; }\n"
                 << "          h2 { margin-top: 0; }\n"
                 << "          h4 { margin-bottom: 5px; margin-top: 15px; color: #555; }\n"
                 << "        </style>\n"
                 << "      </head>\n"
                 << "      <body>\n"
                 << "        <h1>Darwin Report: <xsl:value-of select=\"@file\"/></h1>\n"
                 << "        <div class=\"row\">\n"
                 << "          <div class=\"info\">\n"
                 << "            <h2>Whole Sample</h2>\n"
                 << "            <p><strong>Dimensions:</strong> <xsl:value-of select=\"@dimensions\"/></p>\n"
                 << "            <p><strong>Total Events:</strong> <xsl:value-of select=\"Sample/Summary/@totalEvents\"/></p>\n"
                 << "            <p><strong>Clusters Found:</strong> <xsl:value-of select=\"Sample/Summary/@clustersFound\"/></p>\n"
                 << "          </div>\n"
                 << "          <div class=\"plots\">\n"
                 << "            <xsl:for-each select=\"Sample/Visualizations/Visualization\">\n"
                 << "              <div class=\"plot\"><img src=\"{@src}\"/></div>\n"
                 << "            </xsl:for-each>\n"
                 << "          </div>\n"
                 << "        </div>\n"
                 << "          <xsl:for-each select=\"Clusters/Cluster\">\n"
                 << "          <div class=\"row\">\n"
                 << "            <div class=\"info\">\n"
                 << "              <h2>Cluster <xsl:value-of select=\"@id\"/></h2>\n"
                 << "              <p><strong>Events:</strong> <xsl:value-of select=\"@events\"/> (<xsl:value-of select=\"@percentage\"/>%)</p>\n"
                 << "              <h4>Mean</h4>\n"
                 << "              <table class=\"data-table\">\n"
                 << "                <tr><xsl:for-each select=\"Mean/Value\"><th><xsl:value-of select=\"@label\"/></th></xsl:for-each></tr>\n"
                 << "                <tr><xsl:for-each select=\"Mean/Value\"><td><xsl:value-of select=\".\"/></td></xsl:for-each></tr>\n"
                 << "              </table>\n"
                 << "              <h4>Covariance</h4>\n"
                 << "              <table class=\"data-table\">\n"
                 << "                <xsl:for-each select=\"Covariance/Row\">\n"
                 << "                  <tr><xsl:for-each select=\"Cell\"><td><xsl:value-of select=\".\"/></td></xsl:for-each></tr>\n"
                 << "                </xsl:for-each>\n"
                 << "              </table>\n"
                 << "            </div>\n"
                 << "            <div class=\"plots\">\n"
                 << "              <xsl:for-each select=\"Visualizations/Visualization\">\n"
                 << "                <div class=\"plot\"><img src=\"{@src}\"/></div>\n"
                 << "              </xsl:for-each>\n"
                 << "            </div>\n"
                 << "          </div>\n"
                 << "          </xsl:for-each>\n"
                 << "      </body>\n"
                 << "    </html>\n"
                 << "  </xsl:template>\n"
                 << "</xsl:stylesheet>\n";
    }

    void make_plot(const std::string &path, const std::vector<std::vector<std::vector<double>>> &class_data, const std::vector<std::vector<double>> &quant_data)
    {
        using namespace matplot;

        auto fig = figure(true);
        fig->size(400, 400);
        auto ax = fig->current_axes();
        hold(ax, on);
        // Generate coordinate matrices mapped to the image pixel grid (0 to N-1).
        auto x_range = matplot::linspace(0, quant_data[0].size() - 1, quant_data[0].size());
        auto y_range = matplot::linspace(0, quant_data.size() - 1, quant_data.size());
        auto [X, Y] = matplot::meshgrid(x_range, y_range);
        image(class_data[0], class_data[1], class_data[2]);

        // Overlay contours using the pixel-mapped X and Y ranges
        auto c = contour(ax, X, Y, quant_data);
        c->levels(matplot::iota(0.1, 0.1, 0.9)); // Equivalent to 0.1:0.1:0.9
        c->color("black");
        c->line_width(1.2);
        ax->grid(true);

        // Update limits to fit the data grid
        ax->xlim({0, (double)quant_data[0].size() - 1});
        ax->ylim({0, (double)quant_data.size() - 1});

        ax->xticks({0, .2 * (double)quant_data[0].size(), .4 * (double)quant_data[0].size(), .6 * (double)quant_data[0].size(), .8 * (double)quant_data[0].size(), (double)quant_data[0].size() - 1});
        ax->xticklabels({"0", ".2", ".4", ".6", ".8", "1"});
        
        ax->yticks({0, .2 * (double)quant_data.size(), .4 * (double)quant_data.size(), .6 * (double)quant_data.size(), .8 * (double)quant_data.size(), (double)quant_data.size() - 1});
        ax->yticklabels({"1", ".8", ".6", ".4", ".2", "0"});
        
        fig->save(path);
    }

    void dispatch_plot(const std::string &path, const std::vector<std::vector<std::vector<double>>> &class_data, const std::vector<std::vector<double>> &quant_data) {
        std::cout.flush();
        if (xml_out.is_open()) xml_out.flush();

#ifndef _WIN32
        while (active_pids.size() >= params.max_concurrent) {
            int status;
            pid_t done = waitpid(-1, &status, 0);
            if (done > 0) {
                active_pids.erase(std::remove(active_pids.begin(), active_pids.end(), done), active_pids.end());
            }
        }

        pid_t pid = fork();
        if (pid == 0) {
            make_plot(path, class_data, quant_data);
            _exit(0);
        } else if (pid > 0) {
            active_pids.push_back(pid);
        } else {
            std::cerr << "Warning: fork() failed! Plotting synchronously." << std::endl;
            make_plot(path, class_data, quant_data);
        }
#else
        // Windows MSVC natively does not support fork(). Fallback to synchronous.
        make_plot(path, class_data, quant_data);
#endif
    }

    void wait_all() {
#ifndef _WIN32
        if (!active_pids.empty()) {
            std::cout << "Waiting for " << active_pids.size() << " background plots to finish rendering..." << std::endl;
            for (pid_t pid : active_pids) {
                int status;
                waitpid(pid, &status, 0);
            }
            active_pids.clear();
        }
#endif
    }
};

template <unsigned Dimension>
int do_it(Darwin &darwin)
{
    Params &params = darwin.params;

    std::cout << "Darwin running with" 
              << " filename=" << params.filename << " smooth=" << params.smooth << " threshold=" << params.threshold
              << " dimension=" << params.dimension << " grid=" << params.grid << std::endl;
    std::cout << "   grow=" << (params.grow ? "on" : "off") << " verify=" << (params.verify ? "on" : "off") << " antialias=" << (params.antialias ? "on" : "off") 
              << " verbose=" << (params.verbose ? "on" : "off") << " visual=" << (params.visual ? "on" : "off")
              << " max-clusters=" << params.max_clusters << " min-events=" << params.min_events << std::endl;

    // Load the event data
    typename Weighty<Dimension>::Events events;
    std::string ext = params.ascii ? ".txt" : ".dat";
    std::cout << "Loading events from " << params.filename << ext << "...";
    if (!events.read(params.filename + ext, params.ascii))
    {
        std::cerr << std::endl << "Error: Could not open event file: " << params.filename << ext << std::endl;
        return 1;
    }
    std::cout << " loaded " << events.size() << " events." << std::endl;

    // Riemann object for the whole n-dimensional sample
    Riemann<Dimension> global(params.grid);
    global.visualize = false; // We will save files manually into the dir
    global.antialias = params.antialias;
    global.verify = params.verify;

    // Weighty object for 2-dimensional marginal distributions
    Weighty<2> marginal(params.grid);
    marginal.visualize = false; // We will save files manually into the dir
    marginal.antialias = params.antialias;
    marginal.verify = params.verify;

    // Compute weights for whole sample
    size_t valid_events = 0;
    for (const auto &e : events)
        if (global.event(e)) ++valid_events;

    std::cout << "Initial global analysis found " << valid_events << " valid events...";
    global.Laplace<Dimension>::analyze(params.smooth, params.threshold);
    darwin.num_clusters = global.cluster(params.threshold, params.grow);
    std::cout << " in " << darwin.num_clusters << " clusters." << std::endl;

    // Save the global classification volume once
    global.cluster_id.write(params.out_dir + "/" + params.filename + "_classes.bin");

    darwin.report["total_events"] = events.size();
    darwin.report["num_clusters"] = darwin.num_clusters;

    // Analysis of the whole sample
    std::vector<std::vector<typename Weighty<Dimension>::Event>> cluster_events(params.max_clusters + 1);
    for (const auto &e : events)
    {
        unsigned short c = 0;
        Coordinate<Dimension> coord(global);
        if (global.locate(e, coord))
        {
            size_t idx = (size_t)coord;
            c = static_cast<const Function<Dimension, unsigned short> &>(global.cluster_id)[idx];
        }
        if (c < params.max_clusters + 1)
            cluster_events[c].push_back(e);
        else
            cluster_events[0].push_back(e);
    }
    std::vector<std::string> sample_images;
    for_each_plane<Dimension>([&darwin, &params, &cluster_events, &marginal, &sample_images](unsigned i, unsigned j)
    {
        marginal.reset();
        Coordinate<2> marginal_coord(marginal);
        Weighty<2>::Event marginal_event;
        std::vector<unsigned short> marginal_klass(params.points * params.points, 0);

        for (unsigned short c = 0; c <= darwin.num_clusters; ++c)
        {
            for (auto &e : cluster_events[c])
            {
                marginal_event[0] = e[i];
                marginal_event[1] = e[j];
                marginal.event(marginal_event);
                if (c == 0) continue;

                if (marginal.locate(marginal_event, marginal_coord))
                {
                    size_t idx = (size_t)marginal_coord;
                    unsigned short d = marginal_klass[idx];
                    if (d == 0) marginal_klass[idx] = c;
                    else marginal_klass[idx] = std::min(c, d);
                }
            }
        }
        marginal.prepare(params.smooth);

        if (params.visual) 
        {
            using namespace matplot;

            std::vector<std::vector<std::vector<double>>> class_data(3, std::vector<std::vector<double>>(params.points, std::vector<double>(params.points)));
            std::vector<std::vector<double>> quant_data(params.points, std::vector<double>(params.points));
            for (unsigned y = 0; y < params.points; ++y) {
                for (unsigned x = 0; x < params.points; ++x) {
                    size_t idx = x + y * params.points;
                    for (unsigned i = 0; i < 3; ++i)
                        if (marginal_klass[idx] == 0)
                            class_data[i][y][x] = 255;
                        else
                        {
                            unsigned hue = (255 * marginal_klass[idx] / (std::min(darwin.num_clusters, params.max_clusters)));
                            class_data[i][y][x] = 255 * darwin.colors[hue][i];
                        }
                    quant_data[y][x] = (double)static_cast<const Function<2, float>&>(marginal.quantile)[idx];
                }
            }

            std::string path = params.img_dir + "/sample_" + darwin.labels[i] + "_" + darwin.labels[j] + ".png";
            darwin.dispatch_plot(path, class_data, quant_data);
            sample_images.push_back("images/sample_" + darwin.labels[i] + "_" + darwin.labels[j] + ".png");
    } });
    
    darwin.xml_out << "  <Sample>\n";
    darwin.xml_out << "    <Summary totalEvents=\"" << events.size() << "\" clustersFound=\"" << darwin.num_clusters << "\" />\n";
    darwin.xml_out << "    <Visualizations>\n";
    for (const auto& img : sample_images) {
        darwin.xml_out << "      <Visualization src=\"" << img << "\" />\n";
    }
    darwin.xml_out << "    </Visualizations>\n";
    darwin.xml_out << "  </Sample>\n";
    darwin.xml_out << "  <Clusters>\n";
    
    for (unsigned c = 1; c <= std::min(darwin.num_clusters, params.max_clusters); ++c)
    {
        size_t cluster_event_count = cluster_events[c].size();
        if (cluster_events[c].size() < params.min_events)
        {
            std::cout << "  Skipping small " << cluster_event_count << "events  cluster " << c << "..." << std::endl;
            continue;
        }

        // Riemannian analysis on the cluster subset
        global.reset();

        std::vector<double> mean(Dimension, 0.0);
        std::vector<std::vector<double>> covariance(Dimension, std::vector<double>(Dimension, 0.0));
        for (auto &e : cluster_events[c])
        {
            global.event(e);
            for (unsigned i = 0; i < Dimension; i++)
                mean[i] += e[i];
        }
        for (unsigned i = 0; i < Dimension; i++)
            mean[i] /= cluster_events[c].size();
        for (auto &e : cluster_events[c])
            for (unsigned i = 0; i < Dimension; i++)
                for (unsigned j = 0; j < Dimension; j++)
                    covariance[i][j] += (e[i] - mean[i]) * (e[j] - mean[j]);
        for (unsigned i = 0; i < Dimension; i++)
            for (unsigned j = 0; j < Dimension; j++)
                covariance[i][j] /= (cluster_events[c].size() - 1);

        double pct = 100.0 * cluster_events[c].size() / valid_events;
        std::cout << "Cluster " << c << " contains " << std::fixed << std::setprecision(1) << pct << "% or " << cluster_events[c].size() << " valid events." << std::endl;
        int w = 10;
        std::cout << "+-----+" << std::string(w + 2, '-') << "+" << std::string((w + 1) * Dimension + 1, '-') << "+" << std::endl;
        std::cout << "| Dim | " << std::setw(w) << std::left << "Mean" << " | " << std::setw((w + 1) * Dimension - 1) << std::left << "Covariance" << " |" << std::endl;
        std::cout << "+-----+" << std::string(w + 2, '-') << "+" << std::string((w + 1) * Dimension + 1, '-') << "+" << std::endl;
        for (unsigned i = 0; i < Dimension; i++)
        {
            std::string label_sub = darwin.labels[i].length() > 3 ? darwin.labels[i].substr(0, 3) : darwin.labels[i];
            std::cout << "| " << std::setw(3) << std::right << label_sub << " | " 
                      << std::setw(w) << std::fixed << std::setprecision(4) << mean[i] << " |";
            for (unsigned j = 0; j < Dimension; j++)
                std::cout << " " << std::setw(w) << std::fixed << std::setprecision(4) << covariance[i][j];
            std::cout << " |" << std::endl;
        }
        std::cout << "+-----+" << std::string(w + 2, '-') << "+" << std::string((w + 1) * Dimension + 1, '-') << "+" << std::endl;

        global.Riemann<Dimension>::analyze(params.smooth, params.threshold);

        std::vector<std::string> cluster_images;
        for_each_plane<Dimension>([&darwin, &params, c, &cluster_events, &marginal, &cluster_images](unsigned i, unsigned j)
        {
            Weighty<2>::Event marginal_event;
            Coordinate<2> marginal_coord(marginal);
            std::vector<unsigned short> marginal_klass(params.points * params.points, 0);

            marginal.reset();
            for (auto &e : cluster_events[c])
            {
                marginal_event[0] = e[i];
                marginal_event[1] = e[j];
                marginal.event(marginal_event);
                if (marginal.locate(marginal_event, marginal_coord))
                {
                    size_t idx = (size_t)marginal_coord;
                    marginal_klass[idx] = c;
                }
            }
            marginal.prepare(params.smooth);

            // make cluster marginal plot
            if (params.visual) 
            {
                std::vector<std::vector<std::vector<double>>> class_data(3, std::vector<std::vector<double>>(params.points, std::vector<double>(params.points)));
                std::vector<std::vector<double>> quant_data(params.points, std::vector<double>(params.points));

                for (unsigned y = 0; y < params.points; ++y) {
                    for (unsigned x = 0; x < params.points; ++x) {
                        size_t idx = x + y * params.points;
                        for (unsigned i = 0; i < 3; ++i)
                            if (marginal_klass[idx] == 0)
                                class_data[i][y][x] = 255;
                            else
                            {
                                unsigned hue = (255 * marginal_klass[idx] / (std::min(darwin.num_clusters, params.max_clusters)));
                                class_data[i][y][x] = 255 * darwin.colors[hue][i];
                            }
                        quant_data[y][x] = (double)static_cast<const Function<2, float>&>(marginal.quantile)[idx];
                    }
                }

                // auto class_data = std::make_tuple(std::move(class_r), std::move(class_g), std::move(class_b));
                std::string path = params.img_dir + "/cluster_" + std::to_string(c) + "_" + darwin.labels[i] + "_" + darwin.labels[j] + ".png";
                darwin.dispatch_plot(path, class_data, quant_data);
                cluster_images.push_back("images/cluster_" + std::to_string(c) + "_" + darwin.labels[i] + "_" + darwin.labels[j] + ".png");
        } });

        // Populate JSON report
        json c_info;
        c_info["id"] = c;
        c_info["event_count"] = cluster_event_count;
        c_info["event_percentage"] = pct;
        c_info["differential_error"] = global.differentialError();
        c_info["factor_probability"] = global.factorProbability();
        c_info["mean"] = mean;
        c_info["covariance"] = covariance;
        c_info["images"] = cluster_images;
        darwin.report["clusters"].push_back(c_info);

        // Populate XML data
        darwin.xml_out << "    <Cluster id=\"" << c << "\" events=\"" << cluster_event_count << "\" percentage=\"" << std::fixed << std::setprecision(1) << pct << "\">\n";
        darwin.xml_out << "      <Metrics diffError=\"" << global.differentialError() << "\" ";
        darwin.xml_out << "factorProb=\"" << global.factorProbability() << "\" />\n";
        darwin.xml_out << "      <Mean>\n";
        for (unsigned i = 0; i < Dimension; ++i)
            darwin.xml_out << "        <Value label=\"" << darwin.labels[i] << "\">" << std::fixed << std::setprecision(4) << mean[i] << "</Value>\n";
        darwin.xml_out << "      </Mean>\n";
        darwin.xml_out << "      <Covariance>\n";
        for (unsigned i = 0; i < Dimension; ++i) {
            darwin.xml_out << "        <Row>\n";
            for (unsigned j = 0; j < Dimension; ++j) {
                darwin.xml_out << "          <Cell>" << std::fixed << std::setprecision(4) << covariance[i][j] << "</Cell>\n";
            }
            darwin.xml_out << "        </Row>\n";
        }
        darwin.xml_out << "      </Covariance>\n";
        darwin.xml_out << "      <Visualizations>\n";
        for (const auto& img : cluster_images) {
            darwin.xml_out << "        <Visualization src=\"" << img << "\" />\n";
        }
        darwin.xml_out << "      </Visualizations>\n";
        darwin.xml_out << "    </Cluster>\n";
    }

    darwin.wait_all();

    darwin.xml_out << "  </Clusters>\n";
    darwin.xml_out << "</DarwinReport>\n";
    darwin.xml_out.close();

    std::ofstream json_out(params.out_dir + "/report.json");
    json_out << darwin.report.dump(4);
    json_out.close();

    std::string xml_path = params.out_dir + "/report.xml";
    std::string default_xsl_path = params.out_dir + "/report.xsl";
    std::string default_out_path = params.out_dir + "/index.html";
    
    std::cout << "Generating HTML report..." << std::endl;
    apply_stylesheet(xml_path, default_xsl_path, default_out_path);

    if (!params.style.empty()) {
        fs::path style_path(params.style);
        std::string stem = style_path.stem().string();
        if (fs::path(stem).extension().empty()) {
            stem += ".html";
        }
        std::string custom_out_path = params.out_dir + "/" + stem;
        std::cout << "Generating custom report..." << std::endl;
        apply_stylesheet(xml_path, params.style, custom_out_path);
    }

    std::cout << "Analysis complete. Reports saved and HTML generated." << std::endl;

    return 0;
}

int main(int argc, char *argv[])
{
    Params params;

    // Force Gnuplot to use a non-interactive terminal. 
#ifdef _WIN32
    _putenv_s("GNUTERM", "png");
#else
    setenv("GNUTERM", "png", 1);
#endif

    cxxopts::Options options("Darwin", "Hierarchical Laplacian and Riemannian analysis");

    options.add_options()
    ("f,file", "Input filename", cxxopts::value<std::string>()->default_value("test_data"))
    ("d,dimension", "Dimension", cxxopts::value<unsigned>()->default_value("2"))
    ("g,grid", "Grid resolution", cxxopts::value<unsigned>())
    ("s,smooth", "Smoothing factor", cxxopts::value<float>()->default_value("0.01"))
    ("t,threshold", "Threshold", cxxopts::value<float>()->default_value("0.001"))
    ("visual", "Save visualization files", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("max-clusters", "Maximum number of clusters to report", cxxopts::value<unsigned>()->default_value("10"))
    ("min-events", "Minimum number of events for a cluster to be reported", cxxopts::value<size_t>()->default_value("2")->implicit_value("100"))
    ("labels", "Comma or space separated list of labels", cxxopts::value<std::string>())
    ("v,verbose", "Verbose", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("antialias", "Antialiasing", cxxopts::value<bool>()->default_value("true")->implicit_value("true"))
    ("verify", "Verify consistency", cxxopts::value<bool>()->default_value("true")->implicit_value("true"))
    ("grow", "Grow clusters", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("a,ascii", "Use ASCII data", cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
    ("style", "Path to custom XSLT stylesheet", cxxopts::value<std::string>())
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
    params.smooth = result["smooth"].as<float>();
    params.threshold = result["threshold"].as<float>();
    params.visual = result["visual"].as<bool>();
    params.max_clusters = result["max-clusters"].as<unsigned>();
    params.min_events = result["min-events"].as<size_t>();
    params.antialias = result["antialias"].as<bool>();
    params.verify = result["verify"].as<bool>();
    params.grow = result["grow"].as<bool>();
    params.verbose = result["verbose"].as<bool>();
    params.ascii = result["ascii"].as<bool>();
    if (result.count("labels"))
        params.labels_str = result["labels"].as<std::string>();
    if (result.count("style"))
        params.style = result["style"].as<std::string>();

    params.out_dir = params.filename + ".darwin";
    params.img_dir = params.out_dir + "/images";

    if (!params.style.empty()) {
        std::string xml_path = params.out_dir + "/report.xml";
        if (fs::exists(xml_path)) {
            fs::path style_path(params.style);
            std::string stem = style_path.stem().string();
            if (fs::path(stem).extension().empty()) {
                stem += ".html";
            }
            std::string custom_out_path = params.out_dir + "/" + stem;
            std::cout << "XML already exists. Rendering custom report directly to " << custom_out_path << "..." << std::endl;
            apply_stylesheet(xml_path, params.style, custom_out_path);
            return 0;
        }
    }

    params.max_concurrent = std::max(1u, std::thread::hardware_concurrency());
    if (params.min_events <= 1)
    {
        std::cerr << "Error: min-events must be greater than 1." << std::endl;
        return 1;
    }
    if (result.count("grid"))
        params.grid = result["grid"].as<unsigned>();
    else
        switch (params.dimension)
        {
        case 2:
            params.grid = 256;
            break;
        case 3:
            params.grid = 64;
            break;
        default:
            params.grid = 32;
            break;
        }
    params.points = params.grid + 1;

    Darwin darwin(params);
    darwin.setup();

    switch (params.dimension)
    {
    case 2:
        return do_it<2>(darwin);
    case 3:
        return do_it<3>(darwin);
    case 4:
        return do_it<4>(darwin);
    default:
        std::cerr << "Unsupported dimension." << std::endl;
        return 1;
    }
}
