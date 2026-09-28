#include "Ledger.hpp"

#include <algorithm>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

#include <nlohmann/json.hpp>

namespace fs = std::filesystem;
using json = nlohmann::json;

static std::string escape_html(const std::string &input)
{
    std::string output;
    output.reserve(input.size());
    for (char c : input)
    {
        switch (c)
        {
        case '&':
            output += "&amp;";
            break;
        case '<':
            output += "&lt;";
            break;
        case '>':
            output += "&gt;";
            break;
        case '\"':
            output += "&quot;";
            break;
        default:
            output += c;
            break;
        }
    }
    return output;
}

std::string Ledger::now_iso8601()
{
    const auto now = std::chrono::system_clock::now();
    const std::time_t tt = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &tt);
#else
    localtime_r(&tt, &tm);
#endif
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
    std::ostringstream stamp;
    stamp << std::put_time(&tm, "%Y-%m-%dT%H:%M:%S")
          << '.' << std::setw(3) << std::setfill('0') << static_cast<int>(ms.count());
    return stamp.str();
}

void Ledger::append(const std::string &report_dir, const AnalysisRecord &record)
{
    std::error_code ec;
    fs::create_directories(report_dir, ec);

    std::ofstream out(report_dir + "/analyses.jsonl", std::ios::app);
    if (!out)
        return;

    json j;
    j["timestamp"] = record.timestamp;
    j["workspace"] = record.workspace;
    j["sample"] = record.sample;
    j["population"] = record.population;
    j["analysis"] = record.analysis;
    j["variables"] = record.variables;
    j["stains"] = record.stains;
    j["total_events"] = record.total_events;
    j["valid_clusters"] = record.valid_clusters;
    j["report_file"] = record.report_file;
    out << j.dump() << "\n";
}

std::vector<AnalysisRecord> Ledger::load_all(const std::string &report_dir)
{
    std::vector<AnalysisRecord> records;
    std::ifstream in(report_dir + "/analyses.jsonl");
    if (!in)
        return records;

    std::string line;
    size_t line_number = 0;
    while (std::getline(in, line))
    {
        ++line_number;
        if (line.find_first_not_of(" \t\r\n") == std::string::npos)
            continue;

        try
        {
            json j = json::parse(line);
            AnalysisRecord record;
            record.timestamp = j.value("timestamp", "");
            record.workspace = j.value("workspace", "");
            record.sample = j.value("sample", "");
            record.population = j.value("population", "");
            record.analysis = j.value("analysis", "");
            record.variables = j.value("variables", std::vector<std::string>{});
            record.stains = j.value("stains", std::vector<std::string>{});
            record.total_events = j.value("total_events", static_cast<size_t>(0));
            record.valid_clusters = j.value("valid_clusters", 0u);
            record.report_file = j.value("report_file", "");
            records.push_back(std::move(record));
        }
        catch (const std::exception &e)
        {
            std::cerr << "Ledger: skipping malformed record at " << report_dir << "/analyses.jsonl:"
                       << line_number << " (" << e.what() << ")\n";
        }
    }
    return records;
}

std::string Ledger::stain_set_key(const std::vector<std::string> &stains)
{
    std::vector<std::string> sorted_stains = stains;
    std::sort(sorted_stains.begin(), sorted_stains.end());
    std::string key;
    for (size_t i = 0; i < sorted_stains.size(); ++i)
    {
        if (i)
            key += "+";
        key += sorted_stains[i];
    }
    return key.empty() ? "(none)" : key;
}

static std::string describe(const AnalysisRecord &r)
{
    std::ostringstream ss;
    ss << r.analysis << " (" << r.total_events << " events";
    if (r.valid_clusters)
        ss << ", " << r.valid_clusters << " clusters";
    ss << ") @ " << r.timestamp << " -> " << r.report_file;
    return ss.str();
}

static void print_by_sample(const std::vector<AnalysisRecord> &records)
{
    std::map<std::string, std::map<std::string, std::vector<const AnalysisRecord *>>> grouped;
    for (const auto &r : records)
        grouped[r.sample][r.population].push_back(&r);

    std::cout << "=== Analyses by Sample ===\n";
    for (const auto &[sample, populations] : grouped)
    {
        std::cout << sample << "\n";
        for (const auto &[population, recs] : populations)
        {
            std::cout << "  " << population << "\n";
            for (const auto *r : recs)
                std::cout << "    - " << describe(*r) << "\n";
        }
    }
}

static void print_by_population(const std::vector<AnalysisRecord> &records)
{
    std::map<std::string, std::map<std::string, std::vector<const AnalysisRecord *>>> grouped;
    for (const auto &r : records)
        grouped[r.population][r.sample].push_back(&r);

    std::cout << "=== Analyses by Population ===\n";
    for (const auto &[population, samples] : grouped)
    {
        std::cout << population << "\n";
        for (const auto &[sample, recs] : samples)
        {
            std::cout << "  " << sample << "\n";
            for (const auto *r : recs)
                std::cout << "    - " << describe(*r) << "\n";
        }
    }
}

static void print_by_stain_set(const std::vector<AnalysisRecord> &records)
{
    std::map<std::string, std::vector<const AnalysisRecord *>> grouped;
    for (const auto &r : records)
        grouped[Ledger::stain_set_key(r.stains)].push_back(&r);

    std::cout << "=== Analyses by Stain Set ===\n";
    for (const auto &[stain_set, recs] : grouped)
    {
        std::cout << stain_set << "\n";
        for (const auto *r : recs)
            std::cout << "  " << r->sample << " / " << r->population << " - " << describe(*r) << "\n";
    }
}

void Ledger::print_text_summary(const std::vector<AnalysisRecord> &records, const std::string &mode)
{
    if (records.empty())
    {
        std::cout << "No analyses recorded yet.\n";
        return;
    }

    if (mode == "sample" || mode == "all")
        print_by_sample(records);
    if (mode == "population" || mode == "all")
        print_by_population(records);
    if (mode == "stains" || mode == "all")
        print_by_stain_set(records);
}

static std::string html_table(const std::vector<const AnalysisRecord *> &recs, bool show_sample, bool show_population)
{
    std::ostringstream out;
    out << "<table><tr>";
    if (show_sample)
        out << "<th>Sample</th>";
    if (show_population)
        out << "<th>Population</th>";
    out << "<th>Analysis</th><th>Events</th><th>Clusters</th><th>Timestamp</th><th>Report</th></tr>\n";
    for (const auto *r : recs)
    {
        out << "<tr>";
        if (show_sample)
            out << "<td>" << escape_html(r->sample) << "</td>";
        if (show_population)
            out << "<td>" << escape_html(r->population) << "</td>";
        out << "<td>" << escape_html(r->analysis) << "</td>"
            << "<td>" << r->total_events << "</td>"
            << "<td>" << (r->valid_clusters ? std::to_string(r->valid_clusters) : "-") << "</td>"
            << "<td>" << escape_html(r->timestamp) << "</td>"
            << "<td><a href=\"" << escape_html(r->report_file) << "\">" << escape_html(r->report_file) << "</a></td>"
            << "</tr>\n";
    }
    out << "</table>\n";
    return out.str();
}

static std::string build_sample_panel(const std::vector<AnalysisRecord> &records)
{
    std::map<std::string, std::map<std::string, std::vector<const AnalysisRecord *>>> grouped;
    for (const auto &r : records)
        grouped[r.sample][r.population].push_back(&r);

    std::ostringstream out;
    for (const auto &[sample, populations] : grouped)
    {
        out << "<h3>" << escape_html(sample) << "</h3>\n";
        for (const auto &[population, recs] : populations)
        {
            out << "<h4>" << escape_html(population) << "</h4>\n"
                << html_table(recs, false, false);
        }
    }
    return out.str();
}

static std::string build_population_panel(const std::vector<AnalysisRecord> &records)
{
    std::map<std::string, std::map<std::string, std::vector<const AnalysisRecord *>>> grouped;
    for (const auto &r : records)
        grouped[r.population][r.sample].push_back(&r);

    std::ostringstream out;
    for (const auto &[population, samples] : grouped)
    {
        out << "<h3>" << escape_html(population) << "</h3>\n";
        for (const auto &[sample, recs] : samples)
        {
            out << "<h4>" << escape_html(sample) << "</h4>\n"
                << html_table(recs, false, false);
        }
    }
    return out.str();
}

static std::string build_stain_panel(const std::vector<AnalysisRecord> &records)
{
    std::map<std::string, std::vector<const AnalysisRecord *>> grouped;
    for (const auto &r : records)
        grouped[Ledger::stain_set_key(r.stains)].push_back(&r);

    std::ostringstream out;
    for (const auto &[stain_set, recs] : grouped)
    {
        out << "<h3>" << escape_html(stain_set) << "</h3>\n"
            << html_table(recs, true, true);
    }
    return out.str();
}

std::string Ledger::write_summary_html(const std::string &report_dir, const std::vector<AnalysisRecord> &records)
{
    std::error_code ec;
    fs::create_directories(report_dir, ec);

    const std::string html_path = report_dir + "/summary.html";
    std::ofstream out(html_path);

    out << R"HTML(<!DOCTYPE html>
<html><head>
<title>Analysis History</title>
<style>
    body { font-family: sans-serif; margin: 20px; background: #f4f4f9; color: #333; }
    .container { max-width: 1100px; margin: 0 auto; background: #fff; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
    h1 { border-bottom: 2px solid #eee; padding-bottom: 10px; }
    h3 { margin-top: 20px; color: #1d4ed8; }
    h4 { margin: 10px 0 4px 0; color: #475569; }
    table { border-collapse: collapse; width: 100%; margin-bottom: 10px; }
    th, td { border: 1px solid #ccc; padding: 6px; text-align: center; font-size: 13px; }
    th { background: #f9f9f9; }
    .tabbed input[name="ledger-tab"] { display: none; }
    .tab-bar { display: flex; gap: 4px; border-bottom: 2px solid #e2e8f0; margin-bottom: 16px; }
    .tab-bar label { padding: 8px 16px; cursor: pointer; border-radius: 6px 6px 0 0; background: #f1f5f9; color: #475569; font-weight: 600; }
    .tab-panel { display: none; }
    #tab-sample:checked ~ #panel-sample,
    #tab-population:checked ~ #panel-population,
    #tab-stains:checked ~ #panel-stains { display: block; }
    #tab-sample:checked ~ .tab-bar label[for="tab-sample"],
    #tab-population:checked ~ .tab-bar label[for="tab-population"],
    #tab-stains:checked ~ .tab-bar label[for="tab-stains"] { background: #2563eb; color: #fff; }
</style>
</head><body>
<div class="container">
<h1>Analysis History</h1>
<div class="tabbed">
<input type="radio" name="ledger-tab" id="tab-sample" checked>
<input type="radio" name="ledger-tab" id="tab-population">
<input type="radio" name="ledger-tab" id="tab-stains">
<div class="tab-bar">
<label for="tab-sample">By Sample</label>
<label for="tab-population">By Population</label>
<label for="tab-stains">By Stain Set</label>
</div>
<div class="tab-panel" id="panel-sample">
)HTML"
        << build_sample_panel(records)
        << "</div>\n<div class=\"tab-panel\" id=\"panel-population\">\n"
        << build_population_panel(records)
        << "</div>\n<div class=\"tab-panel\" id=\"panel-stains\">\n"
        << build_stain_panel(records)
        << "</div>\n</div>\n</div>\n</body></html>\n";

    return "summary.html";
}
