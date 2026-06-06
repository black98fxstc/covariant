#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <functional>

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
    };

    std::string id;
    std::string parent_id;
    std::vector<Dimension> dimensions;
    std::vector<std::shared_ptr<Gate>> children;

    Gate(const std::string &id, const std::string &parent_id)
        : id(id), parent_id(parent_id) {}

    virtual ~Gate();

    virtual void apply(DataSet &data);
    virtual void evaluate(const std::vector<std::reference_wrapper<Variable>> &dims, Variable &gate_var) const = 0;
};

class RectangleGate : public Gate
{
public:
    RectangleGate(const std::string &id, const std::string &parent_id)
        : Gate(id, parent_id) {}
    ~RectangleGate() override;

    void evaluate(const std::vector<std::reference_wrapper<Variable>> &dims, Variable &gate_var) const override;
};

class PolygonGate : public Gate
{
public:
    PolygonGate(const std::string &id, const std::string &parent_id)
        : Gate(id, parent_id) {}
    ~PolygonGate() override;

    void evaluate(const std::vector<std::reference_wrapper<Variable>> &dims, Variable &gate_var) const override;
};

class BooleanGate : public Gate
{
public:
    BooleanGate(const std::string &id, const std::string &parent_id)
        : Gate(id, parent_id) {}
    ~BooleanGate() override;

    void evaluate(const std::vector<std::reference_wrapper<Variable>> &dims, Variable &gate_var) const override;
};

class EllipsoidGate : public Gate
{
public:
    EllipsoidGate(const std::string &id, const std::string &parent_id)
        : Gate(id, parent_id) {}
    ~EllipsoidGate() override;

    void evaluate(const std::vector<std::reference_wrapper<Variable>> &dims, Variable &gate_var) const override;
};

class QuadrantGate : public Gate
{
public:
    QuadrantGate(const std::string &id, const std::string &parent_id)
        : Gate(id, parent_id) {}
    ~QuadrantGate() override;

    void evaluate(const std::vector<std::reference_wrapper<Variable>> &dims, Variable &gate_var) const override;
};

void walkGatingTree(const std::shared_ptr<Gate> &node, int depth = 0);