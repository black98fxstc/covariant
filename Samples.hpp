#pragma once

#include <string>
#include <vector>
#include <functional>
#include <variant>
#include <map>

#include "Gating.hpp"
#include "Transforms.hpp"

class DataSet;

enum struct Scale
{
    unknown,
    linear,
    log,
    biexp,
    logicle
};

// Define parameter structs
struct LinearParams
{
    float min_val = 0.0f;
    float max_val = 1024.0f;
};

struct LogParams
{
    float decades = 4.0f;
};

struct BiexpParams
{
    float positive_decades = 4.0f;
    float width = 10.0f;
};

struct LogicleParams
{
    float t = 262144.0f; // Top of scale
    float w = 0.5f;      // Number of decades in the linear region
    float m = 4.5f;      // Number of decades
    float a = 0.0f;      // Additional negative decades
};

// Safe, modern alternative to a raw union
using Parameters = std::variant<std::monostate, LinearParams, LogParams, BiexpParams, LogicleParams>;

class Variable
{
public:
    Scale scale;
    Parameters params;
    std::shared_ptr<Transform> transform = nullptr;
    std::vector<float> unmixing;
    std::shared_ptr<std::vector<float>> data;

    Variable();

    // Construct with specific parameters and automatically deduce the Scale enum
    // Variable(std::string n, Parameters p);

    virtual ~Variable();

    virtual void evaluate();

    double transform_value(double value) const;
    double inverse_transform(double value) const;
};

class Subset
{
public:
    std::shared_ptr<std::vector<bool>> membership;
    std::shared_ptr<Gate> gate;
};

class Classification
{
public:
    std::shared_ptr<std::vector<unsigned short>> classifications;
};

class DataSet
{
    size_t _size = 0;

public:
    std::unordered_map<std::string, Variable> variable;
    std::vector<Variable *> variables;
    std::unordered_map<std::string, Subset> subset;
    std::unordered_map<std::string, Classification> classification;

    DataSet() = default;
    DataSet(const DataSet &) = delete;
    DataSet &operator=(const DataSet &) = delete;
    DataSet(DataSet &&other) noexcept = default;
    DataSet &operator=(DataSet &&other) noexcept = default;

    ~DataSet() = default;

    void for_each_class(const Classification &cls, std::function<void(unsigned short)> func);
    void for_each_member(const Subset &sub, std::function<void(bool)> func);

    // --- Data Input Reading Logic ---
    void setup_float_variables(const std::vector<std::string> &headers);

    bool read_csv(const std::string &filename, char delimiter = ',');
    bool read_xml_xslt(const std::string &xml_file, const std::string &xsl_file, char delimiter = ',');
    bool read_binary(const std::string &filename, const std::vector<std::string> &headers);
    bool read_fcs(const std::string &filename);
    bool read(const std::string &filename, const std::vector<std::string> &headers = {});
    size_t size() const;
};
