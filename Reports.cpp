#include "Reports.hpp"
#include "Leonard.hpp"
#include <fstream>
#include <filesystem>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <sstream>

#include <libxml/parser.h>
#include <libxslt/xslt.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

namespace fs = std::filesystem;

static std::string escape_xml(const std::string& input) {
    std::string output;
    output.reserve(input.size());
    for (char c : input) {
        switch (c) {
            case '&':  output += "&amp;"; break;
            case '<':  output += "&lt;"; break;
            case '>':  output += "&gt;"; break;
            case '\"': output += "&quot;"; break;
            case '\'': output += "&apos;"; break;
            default:   output += c; break;
        }
    }
    return output;
}

void make_marginal_plot(const std::string &path, const std::vector<std::vector<std::vector<double>>> &class_data, const std::vector<std::vector<double>> &quant_data)
{
    using namespace matplot;

    auto backend = std::make_shared<matplot::backend::gnuplot>();
    auto fig = std::make_shared<matplot::figure_type>(true);
    fig->backend(backend);
    fig->add_axes();
    fig->size(400, 400);
    auto ax = fig->current_axes();
    ax->hold(on);
    // Generate coordinate matrices mapped to the image pixel grid (0 to N-1).
    auto x_range = matplot::linspace(0, quant_data[0].size() - 1, quant_data[0].size());
    auto y_range = matplot::linspace(0, quant_data.size() - 1, quant_data.size());
    auto [X, Y] = matplot::meshgrid(x_range, y_range);
    auto img = ax->image(class_data[0], class_data[1], class_data[2]);

    // Overlay contours using the pixel-mapped X and Y ranges
    auto c = ax->contour(X, Y, quant_data);
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

void make_gating_plot(const std::string &path, const std::vector<std::vector<double>> &quant_data, const Measurement H, const Measurement V, const Polygon &polygon)
{
    using namespace matplot;

    auto backend = std::make_shared<matplot::backend::gnuplot>();
    auto fig = std::make_shared<matplot::figure_type>(true);
    fig->backend(backend);
    fig->add_axes();
    fig->size(400, 400);
    auto ax = fig->current_axes();
    ax->hold(on);

    // Generate coordinate matrices mapped to [0, 1]
    auto x_range = matplot::linspace(0.0, 1.0, quant_data[0].size());
    auto y_range = matplot::linspace(0.0, 1.0, quant_data.size());
    auto [X, Y] = matplot::meshgrid(x_range, y_range);

    // Overlay contours
    auto c = ax->contour(X, Y, quant_data);
    c->levels(matplot::iota(0.1, 0.1, 0.9)); // Equivalent to 0.1:0.1:0.9
    c->color("black");
    c->line_width(1.2);

    // Extract and draw the EPP separatrix
    std::vector<double> px, py;
    for (const auto& pt : polygon) {
        px.push_back(pt.x());
        py.push_back(pt.y());
    }
    if (!px.empty()) {
        auto p = ax->plot(px, py);
        p->color("red");
        p->line_width(2.0);
    }

    ax->grid(true);

    // Update limits to exactly fit the [0, 1] scaled domain
    ax->xlim({0.0, 1.0});
    ax->ylim({0.0, 1.0});
    
    fig->save(path);
}

static void apply_stylesheet(const std::string& xml_path, const std::string& xsl_path, const std::string& out_path) {
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

void Reports::update_index(const std::string& report_dir, const std::string& title, const std::vector<ReportLink>& links) {
    fs::create_directories(report_dir);
    std::string index_path = report_dir + "/index.html";
    std::ofstream out(index_path);

    out << "<!DOCTYPE html>\n<html>\n<head>\n"
        << "  <title>" << title << "</title>\n"
        << "  <style>\n"
        << "    body { font-family: sans-serif; margin: 20px; background: #f4f4f9; color: #333; }\n"
        << "    .container { max-width: 800px; margin: 0 auto; background: #fff; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }\n"
        << "    h1 { border-bottom: 2px solid #eee; padding-bottom: 10px; }\n"
        << "    ul.report-list { list-style: none; padding: 0; }\n"
        << "    ul.report-list li { margin: 10px 0; padding: 10px; border: 1px solid #ddd; border-radius: 5px; background: #fafafa; }\n"
        << "    ul.report-list a { text-decoration: none; color: #007bff; font-weight: bold; }\n"
        << "    ul.report-list a:hover { text-decoration: underline; }\n"
        << "    .desc { font-size: 0.9em; color: #666; margin-top: 5px; }\n"
        << "  </style>\n"
        << "</head>\n<body>\n"
        << "  <div class=\"container\">\n"
        << "    <h1>" << title << "</h1>\n"
        << "    <ul class=\"report-list\">\n";

    for (const auto& link : links) {
        out << "      <li>\n"
            << "        <a href=\"" << link.filename << "\">" << link.title << "</a>\n"
            << "        <div class=\"desc\">" << link.description << "</div>\n"
            << "      </li>\n";
    }

    out << "    </ul>\n  </div>\n</body>\n</html>\n";
}

std::string Reports::generate_laplace_report(const std::string& report_dir, const std::string& sample_name, const std::string& pop_name, const Laplace_Results& res, const std::vector<std::string>& selected_vars) {
    fs::create_directories(report_dir);
    std::string stem = "laplace_" + sample_name + "_" + pop_name;
    std::replace(stem.begin(), stem.end(), ' ', '_');
    
    std::string xml_path = report_dir + "/" + stem + ".xml";
    std::string xsl_path = report_dir + "/laplace_report.xsl";
    std::string html_path = report_dir + "/" + stem + ".html";
    
    std::ofstream xml_out(xml_path);
    xml_out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    xml_out << "<?xml-stylesheet type=\"text/xsl\" href=\"laplace_report.xsl\"?>\n";
    xml_out << "<LaplaceReport sample=\"" << escape_xml(sample_name) << "\" population=\"" << escape_xml(pop_name) << "\">\n";
    xml_out << "  <Summary clustersFound=\"" << res.clusters_found << "\"/>\n";
    xml_out << "  <Clusters>\n";
        
    for (size_t i = 1; i <= res.clusters_found; ++i) {
        xml_out << "    <Cluster id=\"" << i << "\" events=\"" << (i < res.cluster_events[i].size() ? res.cluster_events[i].size() : 0) << "\"/>\n";
    }
    
    // darwin.report["total_events"] = events.size();
    // darwin.report["num_clusters"] = darwin.num_clusters;

    // // Analysis of the whole sample
    // std::vector<std::vector<Event<Dimension>>> cluster_events(params.max_clusters + 1);
    // for (const auto &e : events)
    // {
    //     unsigned short c = 0;
    //     Coordinates<Dimension> coord(global);
    //     if (global.locate(e, coord))
    //     {
    //         size_t idx = (size_t)coord;
    //         c = static_cast<const Function<Dimension, unsigned short> &>(global.cluster_id)[idx];
    //     }
    //     if (c < params.max_clusters + 1)
    //         cluster_events[c].push_back(e);
    //     else
    //         cluster_events[0].push_back(e);
    // }
    // std::vector<std::string> sample_images;
    // for_each_plane<Dimension>([&darwin, &params, &cluster_events, &marginal, &sample_images](unsigned i, unsigned j)
    // {
    //     marginal.reset();
    //     Coordinates<2> marginal_coord(marginal);
    //     Event<2> marginal_event;
    //     std::vector<unsigned short> marginal_klass(params.points * params.points, 0);

    //     for (unsigned short c = 0; c <= darwin.num_clusters; ++c)
    //     {
    //         for (auto &e : cluster_events[c])
    //         {
    //             marginal_event[0] = e[i];
    //             marginal_event[1] = e[j];
    //             marginal.event(marginal_event);
    //             if (c == 0) continue;

    //             if (marginal.locate(marginal_event, marginal_coord))
    //             {
    //                 size_t idx = (size_t)marginal_coord;
    //                 unsigned short d = marginal_klass[idx];
    //                 if (d == 0) marginal_klass[idx] = c;
    //                 else marginal_klass[idx] = std::min(c, d);
    //             }
    //         }
    //     }
    //     marginal.prepare(params.smooth);

    //     if (params.visual) 
    //     {
    //         using namespace matplot;

    //         std::vector<std::vector<std::vector<double>>> class_data(3, std::vector<std::vector<double>>(params.points, std::vector<double>(params.points)));
    //         std::vector<std::vector<double>> quant_data(params.points, std::vector<double>(params.points));
    //         for (unsigned y = 0; y < params.points; ++y) {
    //             for (unsigned x = 0; x < params.points; ++x) {
    //                 size_t idx = x + y * params.points;
    //                 for (unsigned i = 0; i < 3; ++i)
    //                     if (marginal_klass[idx] == 0)
    //                         class_data[i][y][x] = 255;
    //                     else
    //                     {
    //                         unsigned hue = (255 * marginal_klass[idx] / (std::min(darwin.num_clusters, params.max_clusters)));
    //                         class_data[i][y][x] = 255 * darwin.colors[hue][i];
    //                     }
    //                 quant_data[y][x] = (double)static_cast<const Function<2, float>&>(marginal.quantile)[idx];
    //             }
    //         }

    //         std::string path = params.img_dir + "/sample_" + darwin.labels[i] + "_" + darwin.labels[j] + ".png";
    //         darwin.dispatch_plot(path, class_data, quant_data);
    //         sample_images.push_back("images/sample_" + darwin.labels[i] + "_" + darwin.labels[j] + ".png");
    // } });
    
    xml_out << "  </Clusters>\n";
    xml_out << "</LaplaceReport>\n";
    xml_out.close();

    if (!fs::exists(xsl_path)) {
        std::ofstream xsl_out(xsl_path);
        xsl_out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                << "<xsl:stylesheet version=\"1.0\" xmlns:xsl=\"http://www.w3.org/1999/XSL/Transform\">\n"
                << "  <xsl:output method=\"html\" indent=\"yes\"/>\n"
                << "  <xsl:template match=\"/LaplaceReport\">\n"
                << "    <html>\n"
                << "      <head>\n"
                << "        <title>Laplace Report: <xsl:value-of select=\"@population\"/></title>\n"
                << "        <style>\n"
                << "          body { font-family: sans-serif; margin: 20px; background: #f4f4f9; color: #333; }\n"
                << "          .container { max-width: 1000px; margin: 0 auto; background: #fff; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }\n"
                << "          table { border-collapse: collapse; width: 100%; margin-top: 20px; }\n"
                << "          th, td { border: 1px solid #ccc; padding: 8px; text-align: center; }\n"
                << "          th { background: #f9f9f9; }\n"
                << "        </style>\n"
                << "      </head>\n"
                << "      <body>\n"
                << "        <div class=\"container\">\n"
                << "          <h2>Laplace Analysis: <xsl:value-of select=\"@sample\"/> / <xsl:value-of select=\"@population\"/></h2>\n"
                << "          <p><strong>Clusters Found:</strong> <xsl:value-of select=\"Summary/@clustersFound\"/></p>\n"
                << "          <table>\n"
                << "            <tr><th>Cluster ID</th><th>Events</th></tr>\n"
                << "            <xsl:for-each select=\"Clusters/Cluster\">\n"
                << "              <tr><td><xsl:value-of select=\"@id\"/></td><td><xsl:value-of select=\"@events\"/></td></tr>\n"
                << "            </xsl:for-each>\n"
                << "          </table>\n"
                << "        </div>\n"
                << "      </body>\n"
                << "    </html>\n"
                << "  </xsl:template>\n"
                << "</xsl:stylesheet>\n";
        xsl_out.close();
    }

    apply_stylesheet(xml_path, xsl_path, html_path);

    return stem + ".html";
}

std::string Reports::generate_epp_report(const std::string& report_dir, const std::string& sample_name, const std::string& pop_name, const Pursuit_Results& res, const std::vector<std::string>& selected_vars) {
    fs::create_directories(report_dir);
    std::string stem = "epp_" + sample_name + "_" + pop_name;
    std::replace(stem.begin(), stem.end(), ' ', '_');
    
    std::string xml_path = report_dir + "/" + stem + ".xml";
    std::string xsl_path = report_dir + "/epp_report.xsl";
    std::string html_path = report_dir + "/" + stem + ".html";

    std::ofstream xml_out(xml_path);
    xml_out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    xml_out << "<?xml-stylesheet type=\"text/xsl\" href=\"epp_report.xsl\"?>\n";
    xml_out << "<EPPReport sample=\"" << escape_xml(sample_name) << "\" population=\"" << escape_xml(pop_name) << "\">\n";
    xml_out << "  <AllEvents events=\"" << res.event_count << "\">\n";
    if (!res.means.empty()) {
        xml_out << "    <Means>\n";
        for (size_t i = 0; i < selected_vars.size() && i < res.means.size(); ++i) {
            std::ostringstream ss, ss_pct;
            ss << std::fixed << std::setprecision(4) << res.means[i];
            double pct = std::max(0.0, std::min(100.0, res.means[i] * 100.0));
            ss_pct << std::fixed << std::setprecision(1) << pct;
            xml_out << "      <Variable name=\"" << escape_xml(selected_vars[i]) << "\" mean=\"" << ss.str() << "\" pct=\"" << ss_pct.str() << "\"/>\n";
        }
        xml_out << "    </Means>\n";
    }
    for (const auto& child : res.children) {
        render_epp_node(xml_out, child, selected_vars, 2);
    }
    xml_out << "  </AllEvents>\n";
    xml_out << "</EPPReport>\n";
    xml_out.close();

    std::ofstream xsl_out(xsl_path);
    xsl_out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
            << "<xsl:stylesheet version=\"1.0\" xmlns:xsl=\"http://www.w3.org/1999/XSL/Transform\">\n"
            << "  <xsl:output method=\"html\" indent=\"yes\"/>\n"
            << "  <xsl:template match=\"/EPPReport\">\n"
            << "    <html>\n"
            << "      <head>\n"
            << "        <title>EPP Report: <xsl:value-of select=\"@population\"/></title>\n"
            << "        <style>\n"
            << "          body { font-family: -apple-system, BlinkMacSystemFont, \"Segoe UI\", Roboto, Helvetica, Arial, sans-serif; margin: 20px; background: #f8fafc; color: #1e293b; }\n"
            << "          .container { width: 100%; min-width: 800px; margin: 0 auto; background: #fff; padding: 24px; border-radius: 8px; overflow-x: auto; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1), 0 2px 4px -2px rgba(0,0,0,0.1); box-sizing: border-box; }\n"
            << "          h2 { margin-top: 0; color: #0f172a; font-size: 20px; border-bottom: 2px solid #e2e8f0; padding-bottom: 10px; }\n"
            << "          .tree ul { padding-top: 20px; position: relative; display: flex; justify-content: center; padding-left: 0; }\n"
            << "          .tree li { float: left; text-align: center; list-style-type: none; position: relative; padding: 20px 8px 0 8px; }\n"
            << "          .tree li::before, .tree li::after { content: ''; position: absolute; top: 0; right: 50%; border-top: 2px solid #cbd5e1; width: 50%; height: 20px; }\n"
            << "          .tree li::after { right: auto; left: 50%; border-left: 2px solid #cbd5e1; }\n"
            << "          .tree li:only-child::after, .tree li:only-child::before { display: none; }\n"
            << "          .tree li:only-child { padding-top: 0; }\n"
            << "          .tree li:first-child::before, .tree li:last-child::after { border: 0 none; }\n"
            << "          .tree li:last-child::before { border-right: 2px solid #cbd5e1; border-radius: 0 6px 0 0; }\n"
            << "          .tree li:first-child::after { border-radius: 6px 0 0 0; }\n"
            << "          .tree ul ul::before { content: ''; position: absolute; top: 0; left: 50%; border-left: 2px solid #cbd5e1; width: 0; height: 20px; }\n"
            << "          .tree li div.node { border: 1px solid #cbd5e1; padding: 12px; text-decoration: none; color: #1e293b; font-size: 12px; display: inline-block; border-radius: 6px; background-color: #ffffff; box-shadow: 0 2px 5px rgba(0,0,0,0.06); min-width: 170px; max-width: 220px; vertical-align: top; box-sizing: border-box; }\n"
            << "          .tree li div.node-root { background-color: #f0f7ff; border: 2px solid #3b82f6; }\n"
            << "          .tree li div.node h4 { margin: 2px 0 6px 0; color: #1d4ed8; font-size: 13px; font-weight: 700; word-break: break-word; }\n"
            << "          .tree li div.node p { margin: 3px 0; color: #475569; }\n"
            << "          .tree li div.node img { max-width: 160px; width: 100%; height: auto; display: block; margin: 6px auto; border: 1px solid #e2e8f0; border-radius: 4px; background: #fff; }\n"
            << "          .badge { display: inline-block; padding: 2px 6px; font-size: 10px; font-weight: 700; border-radius: 4px; margin-bottom: 6px; margin-right: 3px; letter-spacing: 0.3px; }\n"
            << "          .badge-in { background-color: #dcfce7; color: #15803d; border: 1px solid #86efac; }\n"
            << "          .badge-out { background-color: #fee2e2; color: #b91c1c; border: 1px solid #fca5a5; }\n"
            << "          .badge-leaf { background-color: #fef3c7; color: #b45309; border: 1px solid #fcd34d; }\n"
            << "          .plot-box { text-align: center; font-size: 10px; color: #64748b; margin: 6px 0; }\n"
            << "          .means-panel { margin-top: 10px; padding-top: 8px; border-top: 1px solid #e2e8f0; text-align: left; }\n"
            << "          .means-title { font-size: 10px; font-weight: 700; color: #64748b; text-transform: uppercase; letter-spacing: 0.5px; margin-bottom: 6px; text-align: center; }\n"
            << "          .means-bars { display: flex; flex-direction: column; gap: 4px; }\n"
            << "          .bar-row { display: flex; align-items: center; font-size: 10px; line-height: 1.2; }\n"
            << "          .bar-label { width: 58px; font-weight: 600; color: #334155; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; padding-right: 4px; }\n"
            << "          .bar-track { flex: 1; height: 8px; background-color: #e2e8f0; border-radius: 4px; position: relative; overflow: hidden; box-shadow: inset 0 1px 2px rgba(0,0,0,0.08); }\n"
            << "          .bar-fill { height: 100%; background: linear-gradient(90deg, #3b82f6 0%, #2563eb 100%); border-radius: 4px; }\n"
            << "          .bar-val { width: 38px; font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace; font-size: 9.5px; color: #475569; text-align: right; padding-left: 4px; }\n"
            << "        </style>\n"
            << "      </head>\n"
            << "      <body>\n"
            << "        <div class=\"container\">\n"
            << "          <h2>EPP Analysis: <xsl:value-of select=\"@sample\"/> / <xsl:value-of select=\"@population\"/></h2>\n"
            << "          <div class=\"tree\">\n"
            << "            <ul>\n"
            << "              <xsl:apply-templates select=\"AllEvents\"/>\n"
            << "            </ul>\n"
            << "          </div>\n"
            << "        </div>\n"
            << "      </body>\n"
            << "    </html>\n"
            << "  </xsl:template>\n"
            << "  <xsl:template match=\"AllEvents\">\n"
            << "    <li>\n"
            << "      <div class=\"node node-root\">\n"
            << "        <h4>All Events</h4>\n"
            << "        <p><strong>Total Events:</strong> <xsl:value-of select=\"@events\"/></p>\n"
            << "        <xsl:if test=\"not(Node) and Means/Variable\">\n"
            << "          <div class=\"means-panel\">\n"
            << "            <div class=\"means-title\">Expression Levels</div>\n"
            << "            <div class=\"means-bars\">\n"
            << "              <xsl:for-each select=\"Means/Variable\">\n"
            << "                <div class=\"bar-row\">\n"
            << "                  <span class=\"bar-label\" title=\"{@name}\"><xsl:value-of select=\"@name\"/></span>\n"
            << "                  <div class=\"bar-track\" title=\"Mean: {@mean} ({@pct}%)\">\n"
            << "                    <div class=\"bar-fill\" style=\"width: {@pct}%;\"></div>\n"
            << "                  </div>\n"
            << "                  <span class=\"bar-val\"><xsl:value-of select=\"format-number(@mean, '#.##')\"/></span>\n"
            << "                </div>\n"
            << "              </xsl:for-each>\n"
            << "            </div>\n"
            << "          </div>\n"
            << "        </xsl:if>\n"
            << "      </div>\n"
            << "      <xsl:if test=\"Node\">\n"
            << "        <ul>\n"
            << "          <xsl:apply-templates select=\"Node\"/>\n"
            << "        </ul>\n"
            << "      </xsl:if>\n"
            << "    </li>\n"
            << "  </xsl:template>\n"
            << "  <xsl:template match=\"Node\">\n"
            << "    <li>\n"
            << "      <div class=\"node\">\n"
            << "        <xsl:if test=\"@branch = 'in'\">\n"
            << "          <span class=\"badge badge-in\">IN</span>\n"
            << "        </xsl:if>\n"
            << "        <xsl:if test=\"@branch = 'out'\">\n"
            << "          <span class=\"badge badge-out\">OUT</span>\n"
            << "        </xsl:if>\n"
            << "        <xsl:if test=\"@isLeaf = 'true' or not(Node)\">\n"
            << "          <span class=\"badge badge-leaf\">LEAF</span>\n"
            << "        </xsl:if>\n"
            << "        <h4>\n"
            << "          <xsl:choose>\n"
            << "            <xsl:when test=\"@gateX and @gateY\">\n"
            << "              <xsl:value-of select=\"@gateX\"/> vs <xsl:value-of select=\"@gateY\"/>\n"
            << "            </xsl:when>\n"
            << "            <xsl:otherwise>\n"
            << "              Subset <xsl:value-of select=\"@id\"/>\n"
            << "            </xsl:otherwise>\n"
            << "          </xsl:choose>\n"
            << "        </h4>\n"
            << "        <p><strong>Events:</strong> <xsl:value-of select=\"@events\"/></p>\n"
            << "        <xsl:if test=\"@image\">\n"
            << "          <div class=\"plot-box\">\n"
            << "            <a href=\"{@image}\" target=\"_blank\">\n"
            << "              <img src=\"{@image}\" alt=\"Gate\" title=\"Click to enlarge\"/>\n"
            << "            </a>\n"
            << "          </div>\n"
            << "        </xsl:if>\n"
            << "        <xsl:if test=\"(@isLeaf = 'true' or not(Node)) and Means/Variable\">\n"
            << "          <div class=\"means-panel\">\n"
            << "            <div class=\"means-title\">Expression Levels</div>\n"
            << "            <div class=\"means-bars\">\n"
            << "              <xsl:for-each select=\"Means/Variable\">\n"
            << "                <div class=\"bar-row\">\n"
            << "                  <span class=\"bar-label\" title=\"{@name}\"><xsl:value-of select=\"@name\"/></span>\n"
            << "                  <div class=\"bar-track\" title=\"Mean: {@mean} ({@pct}%)\">\n"
            << "                    <div class=\"bar-fill\" style=\"width: {@pct}%;\"></div>\n"
            << "                  </div>\n"
            << "                  <span class=\"bar-val\"><xsl:value-of select=\"format-number(@mean, '#.##')\"/></span>\n"
            << "                </div>\n"
            << "              </xsl:for-each>\n"
            << "            </div>\n"
            << "          </div>\n"
            << "        </xsl:if>\n"
            << "      </div>\n"
            << "      <xsl:if test=\"Node\">\n"
            << "        <ul>\n"
            << "          <xsl:apply-templates select=\"Node\"/>\n"
            << "        </ul>\n"
            << "      </xsl:if>\n"
            << "    </li>\n"
            << "  </xsl:template>\n"
            << "</xsl:stylesheet>\n";
    xsl_out.close();

    apply_stylesheet(xml_path, xsl_path, html_path);

    return stem + ".html";
}

void Reports::render_epp_node(std::ostream& out, const Pursuit_Results& node, const std::vector<std::string>& selected_vars, int depth) {
    std::string indent(depth * 2, ' ');
    std::string x_name = (node.has_gate && node.gate_x < selected_vars.size()) ? selected_vars[node.gate_x] : "X";
    std::string y_name = (node.has_gate && node.gate_y < selected_vars.size()) ? selected_vars[node.gate_y] : "Y";

    out << indent << "<Node id=\"" << escape_xml(node.node_id) 
        << "\" branch=\"" << escape_xml(node.branch) 
        << "\" events=\"" << node.event_count 
        << "\" isLeaf=\"" << (node.is_leaf ? "true" : "false") << "\"";

    if (node.has_gate) {
        out << " gateX=\"" << escape_xml(x_name) << "\" gateY=\"" << escape_xml(y_name) << "\"";
    }
    if (!node.polygon_image.empty()) {
        out << " image=\"" << escape_xml(node.polygon_image) << "\"";
    }
    out << ">\n";

    if (!node.means.empty()) {
        out << indent << "  <Means>\n";
        for (size_t i = 0; i < selected_vars.size() && i < node.means.size(); ++i) {
            std::ostringstream ss, ss_pct;
            ss << std::fixed << std::setprecision(4) << node.means[i];
            double pct = std::max(0.0, std::min(100.0, node.means[i] * 100.0));
            ss_pct << std::fixed << std::setprecision(1) << pct;
            out << indent << "    <Variable name=\"" << escape_xml(selected_vars[i]) << "\" mean=\"" << ss.str() << "\" pct=\"" << ss_pct.str() << "\"/>\n";
        }
        out << indent << "  </Means>\n";
    }

    if (!node.children.empty()) {
        for (const auto& child : node.children) {
            render_epp_node(out, child, selected_vars, depth + 1);
        }
    }

    out << indent << "</Node>\n";
}