#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <limits>

#include "Transforms.hpp"

class DataSet;
class Variable;

// ----------------------------------------------------------------------------
// Gating-ML Gate Base Class and Stubs
// ----------------------------------------------------------------------------

class Gate
{
public:
    struct Dimension
    {
        std::string name;
        std::string scale;
        std::string compensation;
        std::shared_ptr<Transform> transform;
        double min_val = -std::numeric_limits<double>::infinity();
        double max_val = std::numeric_limits<double>::infinity();
    };

    std::string id;
    std::string parent_id;
    std::vector<Dimension> dimensions;
    std::vector<std::shared_ptr<std::vector<float>>> data;
    std::shared_ptr<std::vector<bool>> membership;
    std::vector<std::shared_ptr<Gate>> children;

    Gate(const std::string &id, const std::string &parent_id)
        : id(id), parent_id(parent_id) {}

    virtual ~Gate();

    virtual void apply(std::vector<bool> &membership);
    virtual void evaluate(std::vector<bool> &membership) const = 0;
};

class RectangleGate : public Gate
{
public:
    // min and max bounds are stored per dimension inside Gate::dimensions

    RectangleGate(const std::string &id, const std::string &parent_id)
        : Gate(id, parent_id) {}
    ~RectangleGate() override;

    void evaluate(std::vector<bool> &membership) const override;
};

class PolygonGate : public Gate
{
public:
    std::vector<std::vector<double>> vertices;

    PolygonGate(const std::string &id, const std::string &parent_id)
        : Gate(id, parent_id) {}
    ~PolygonGate() override;

    void evaluate(std::vector<bool> &membership) const override;
};

class BooleanGate : public Gate
{
public:
    BooleanGate(const std::string &id, const std::string &parent_id)
        : Gate(id, parent_id) {}
    ~BooleanGate() override;

    void evaluate(std::vector<bool> &membership) const override;
};

class EllipsoidGate : public Gate
{
public:
    std::vector<double> mean;
    std::vector<std::vector<double>> covariance_matrix;

    EllipsoidGate(const std::string &id, const std::string &parent_id)
        : Gate(id, parent_id) {}
    ~EllipsoidGate() override;

    void evaluate(std::vector<bool> &membership) const override;
};

class QuadrantGate : public Gate
{
public:
    QuadrantGate(const std::string &id, const std::string &parent_id)
        : Gate(id, parent_id) {}
    ~QuadrantGate() override;

    void evaluate(std::vector<bool> &membership) const override;
};

void walkGatingTree(const std::shared_ptr<Gate> &node, int depth = 0);