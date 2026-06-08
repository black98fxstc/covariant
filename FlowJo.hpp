#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>

#include "Gating.hpp"
#include "Transforms.hpp"

extern std::vector<std::string> analysis_choices;

struct SpilloverMatrix
{
    std::string id;
    std::string name;
    std::string prefix;
    std::string suffix;
    std::vector<std::string> parameters;
    std::unordered_map<std::string, std::string> comp_infix_map;
    std::vector<std::vector<double>> matrix;
};

struct SampleData
{
    std::string name;
    bool selected = false;
    bool enabled = true;
    std::vector<std::string> variables;
    std::vector<std::string> populations;
    std::vector<std::string> stains;
    std::vector<std::shared_ptr<Gate>> gates;
    std::unordered_map<std::string, std::shared_ptr<Transform>> transforms;
    SpilloverMatrix spillover_matrix;
    std::unordered_map<std::string, std::string> pop_to_gate_id;
    std::unordered_map<std::string, std::string> gate_id_to_pop;
};

struct Workspace
{
    std::string filename;
    std::vector<SampleData> samples;
    std::vector<std::string> all_populations;
    std::vector<std::string> all_variables;
};

struct SelectionState
{
    bool cancelled = false;
    std::vector<SampleData*> samples;
    std::vector<std::string> variables;
    std::vector<std::string> populations;
    int analysis_choice = 0;
    std::string smoothing_str;
    std::string threshold_str;
    std::string max_clusters_str;
    std::string min_events_str;
    float smoothing = 0.01f;
    float threshold = 0.001f;
    unsigned max_clusters = 50;
    size_t min_events = 100;
    unsigned grid_size = 256;
};

std::string find_workspace(int argc, char *argv[]);
Workspace parse_workspace(const std::string &filename);
SelectionState build_ftxui_interface(Workspace &ws);