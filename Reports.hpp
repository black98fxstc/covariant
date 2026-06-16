#pragma once

#include <string>
#include <vector>
#include "Leonard.hpp"

class Reports {
public:
    struct ReportLink {
        std::string title;
        std::string filename;
        std::string description;
    };

    static void update_index(const std::string& report_dir, const std::string& title, const std::vector<ReportLink>& links);

    static std::string generate_laplace_report(const std::string& report_dir, const std::string& sample_name, const std::string& pop_name, const Laplace_Results& res, const std::vector<std::string>& selected_vars);

    static std::string generate_epp_report(const std::string& report_dir, const std::string& sample_name, const std::string& pop_name, const Pursuit_Results& res, const std::vector<std::string>& selected_vars);

private:
    static void render_epp_node(std::ostream& out, const Pursuit_Results& node, const std::vector<std::string>& selected_vars, int depth, int id);
};