#pragma once

#include <string>
#include <vector>
#include <functional>
#include <variant>

enum struct Scale { unknown, linear, log, biexp, logicle };

// Define parameter structs
struct LinearParams {
    float min_val = 0.0f;
    float max_val = 1024.0f;
};

struct LogParams {
    float decades = 4.0f;
};

struct BiexpParams {
    float positive_decades = 4.0f;
    float width = 10.0f;
};

struct LogicleParams {
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
    std::string name;
    Scale scale;
    Parameters params;

    // Using a variant of vectors (Structure of Arrays) is much more efficient and type-safe.
    // It avoids padding overhead per element and provides great cache locality.
    std::variant<std::vector<float>, std::vector<unsigned short>, std::vector<bool>> data;

    Variable(std::string n);

    // Construct with specific parameters and automatically deduce the Scale enum
    Variable(std::string n, Parameters p);

    // Template helper to safely get data of a specific type
    template <typename T>
    std::vector<T>& get_data()
    {
        return std::get<std::vector<T>>(data);
    }

    template <typename T>
    const std::vector<T>& get_data() const
    {
        return std::get<std::vector<T>>(data);
    }
};

// Base interface for polymorphism
class IDataSet
{
public:
    virtual ~IDataSet() = default;
    virtual Variable& operator[](size_t i) = 0;
    virtual const Variable& operator[](size_t i) const = 0;
    virtual size_t size() const = 0;
};

class DataSet : public IDataSet
{
    std::vector<Variable> variables;

public:
    void add_variable(const std::string& name);
    void add_classification(const std::string& name);
    void add_data(const std::vector<float>& datum);

    void for_each_class(const Variable& classification, std::function<void(unsigned short)> func);
    void for_each_member(const Variable& membership, std::function<void(bool)> func);

    // --- Data Input Reading Logic ---
    void setup_float_variables(const std::vector<std::string>& headers);

    bool read_csv(const std::string& filename, char delimiter = ',');
    bool read_xml_xslt(const std::string& xml_file, const std::string& xsl_file, char delimiter = ',');
    bool read_binary(const std::string& filename, const std::vector<std::string>& headers);
    bool read_fcs(const std::string& filename);
    bool read(const std::string& filename, const std::vector<std::string>& headers = {});

    Variable& operator[](size_t i) override;
    const Variable& operator[](size_t i) const override;
    size_t size() const override;
};

class Projection : public IDataSet
{
    std::vector<size_t> variable_indices;
    DataSet& base_dataset;

public:
    Projection(std::vector<size_t> indices, DataSet& base);

    Variable& operator[](size_t i) override;
    const Variable& operator[](size_t i) const override;
    size_t size() const override;
};