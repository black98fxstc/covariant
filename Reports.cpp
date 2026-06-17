#include "Reports.hpp"
#include <fstream>
#include <filesystem>
#include <iostream>
#include <algorithm>

#include <libxml/parser.h>
#include <libxslt/xslt.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

namespace fs = std::filesystem;

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
    xml_out << "<LaplaceReport sample=\"" << sample_name << "\" population=\"" << pop_name << "\">\n";
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
    xml_out << "<EPPReport sample=\"" << sample_name << "\" population=\"" << pop_name << "\">\n";
    render_epp_node(xml_out, res, selected_vars, 0, 0);
    xml_out << "</EPPReport>\n";
    xml_out.close();

    if (!fs::exists(xsl_path)) {
        std::ofstream xsl_out(xsl_path);
        xsl_out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                << "<xsl:stylesheet version=\"1.0\" xmlns:xsl=\"http://www.w3.org/1999/XSL/Transform\">\n"
                << "  <xsl:output method=\"html\" indent=\"yes\"/>\n"
                << "  <xsl:template match=\"/EPPReport\">\n"
                << "    <html>\n"
                << "      <head>\n"
                << "        <title>EPP Report: <xsl:value-of select=\"@population\"/></title>\n"
                << "        <style>\n"
                << "          body { font-family: sans-serif; margin: 20px; background: #f4f4f9; color: #333; }\n"
                << "          .container { width: 100%; min-width: 800px; margin: 0 auto; background: #fff; padding: 20px; border-radius: 8px; overflow-x: auto; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }\n"
                << "          .tree ul { padding-top: 20px; position: relative; display: flex; justify-content: center; padding-left: 0; }\n"
                << "          .tree li { float: left; text-align: center; list-style-type: none; position: relative; padding: 20px 5px 0 5px; }\n"
                << "          .tree li::before, .tree li::after{ content: ''; position: absolute; top: 0; right: 50%; border-top: 2px solid #ccc; width: 50%; height: 20px; }\n"
                << "          .tree li::after{ right: auto; left: 50%; border-left: 2px solid #ccc; }\n"
                << "          .tree li:only-child::after, .tree li:only-child::before { display: none; }\n"
                << "          .tree li:only-child{ padding-top: 0; }\n"
                << "          .tree li:first-child::before, .tree li:last-child::after{ border: 0 none; }\n"
                << "          .tree li:last-child::before{ border-right: 2px solid #ccc; border-radius: 0 5px 0 0; }\n"
                << "          .tree li:first-child::after{ border-radius: 5px 0 0 0; }\n"
                << "          .tree ul ul::before{ content: ''; position: absolute; top: 0; left: 50%; border-left: 2px solid #ccc; width: 0; height: 20px; }\n"
                << "          .tree li div.node { border: 1px solid #ccc; padding: 10px; text-decoration: none; color: #333; font-size: 12px; display: inline-block; border-radius: 5px; background-color: white; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }\n"
                << "          .tree li div.node h4 { margin: 0 0 5px 0; color: #007bff; }\n"
                << "          .tree li div.node p { margin: 2px 0; }\n"
                << "          .tree li div.node img { max-width: 150px; display: block; margin: 5px auto; border: 1px solid #ddd; }\n"
                << "        </style>\n"
                << "      </head>\n"
                << "      <body>\n"
                << "        <div class=\"container\">\n"
                << "          <h2>EPP Analysis: <xsl:value-of select=\"@sample\"/> / <xsl:value-of select=\"@population\"/></h2>\n"
                << "          <div class=\"tree\">\n"
                << "            <ul>\n"
                << "              <xsl:apply-templates select=\"Node\"/>\n"
                << "            </ul>\n"
                << "          </div>\n"
                << "        </div>\n"
                << "      </body>\n"
                << "    </html>\n"
                << "  </xsl:template>\n"
                << "  <xsl:template match=\"Node\">\n"
                << "    <li>\n"
                << "      <div class=\"node\">\n"
                << "        <xsl:choose>\n"
                << "          <xsl:when test=\"Split\">\n"
                << "            <h4><xsl:value-of select=\"Split/@x\"/> vs <xsl:value-of select=\"Split/@y\"/></h4>\n"
                << "            <p>Score: <xsl:value-of select=\"Split/@score\"/></p>\n"
                << "            <p>In: <xsl:value-of select=\"Split/@in\"/> | Out: <xsl:value-of select=\"Split/@out\"/></p>\n"
                << "          </xsl:when>\n"
                << "          <xsl:otherwise>\n"
                << "            <h4>Leaf Node</h4>\n"
                << "            <p>Terminated</p>\n"
                << "          </xsl:otherwise>\n"
                << "        </xsl:choose>\n"
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
    }

    apply_stylesheet(xml_path, xsl_path, html_path);

    return stem + ".html";
}

void Reports::render_epp_node(std::ostream& out, const Pursuit_Results& node, const std::vector<std::string>& selected_vars, int depth, int id) {
    std::string indent(depth * 2 + 2, ' ');
    out << indent << "<Node id=\"" << id << "\">\n";

    if (node.best_split && node.best_split->outcome == Projection_Results::Status::EPP_success) {
        auto& split = *node.best_split;
        std::string x_name = split.X < selected_vars.size() ? selected_vars[split.X] : "X";
        std::string y_name = split.Y < selected_vars.size() ? selected_vars[split.Y] : "Y";
        out << indent << "  <Split x=\"" << x_name << "\" y=\"" << y_name 
            << "\" score=\"" << split.score 
            << "\" in=\"" << split.in.count 
            << "\" out=\"" << split.out.count << "\">\n";
        out << indent << "    <Polygon>\n";
        for (const auto& point : split.separatrix) {
            out << indent << "      <Vertex x=\"" << point.x() << "\" y=\"" << point.y() << "\"/>\n";
        }
        out << indent << "    </Polygon>\n";
        out << indent << "  </Split>\n";
    }
    
    if (!node.children.empty()) {
        int child_id = 0;
        for (const auto& child : node.children) {
            render_epp_node(out, child, selected_vars, depth + 1, (id * 10) + child_id++);
        }
    }
    
    out << indent << "</Node>\n";
}