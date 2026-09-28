#pragma once

#include <cstddef>
#include <string>
#include <vector>

// One completed sample x population x analysis-method run, appended to
// <workspace>.len/analyses.jsonl so history survives across invocations.
struct AnalysisRecord
{
    std::string timestamp;
    std::string workspace;
    std::string sample;
    std::string population;
    std::string analysis;
    std::vector<std::string> variables;
    std::vector<std::string> stains;
    size_t total_events = 0;
    unsigned valid_clusters = 0;
    std::string report_file;
};

class Ledger
{
public:
    static std::string now_iso8601();

    static void append(const std::string &report_dir, const AnalysisRecord &record);
    static std::vector<AnalysisRecord> load_all(const std::string &report_dir);

    // Grouping key for "stain set": sorted stain names joined with '+'.
    static std::string stain_set_key(const std::vector<std::string> &stains);

    // mode is one of "sample", "population", "stains", or "all".
    static void print_text_summary(const std::vector<AnalysisRecord> &records, const std::string &mode);

    // Writes report_dir/summary.html (tabbed: by sample, by population, by stain set) and returns its filename.
    static std::string write_summary_html(const std::string &report_dir, const std::vector<AnalysisRecord> &records);
};
