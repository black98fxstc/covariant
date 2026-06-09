#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>
#include <algorithm>
#include <cassert>

#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <libxml/tree.h>

#include "Gating.hpp"
#include "Samples.hpp"

// ----------------------------------------------------------------------------
// Parsing and Tree Construction
// ----------------------------------------------------------------------------

Gate::~Gate() = default;

void Gate::apply(std::vector<bool> &membership)
{
    this->evaluate(membership);
}

RectangleGate::~RectangleGate() = default;

void RectangleGate::evaluate(std::vector<bool> &membership)
{
    std::cout << "Applying RectangleGate: " << name ;
    size_t count = 0;
    for (size_t i = 0; i < dimensions.size(); ++i)
    {
        float *bare_data = data[i].get()->data();
        float lower = dimensions[i].transform->scale(dimensions[i].min_val);
        for (size_t j = 0; j < data[i]->size(); j++)
            if (bare_data[j] < lower)
                membership[j] = false;
        float upper = dimensions[i].transform->scale(dimensions[i].max_val);
        for (size_t j = 0; j < data[i]->size(); j++)
            if (bare_data[j] >= upper)
                membership[j] = false;
    }
    for (size_t j = 0; j < data[0]->size(); j++)
        if (membership[j])
            ++count;
    std::cout << " found " << count << " events." << std::endl;
}

PolygonGate::~PolygonGate() = default;

void PolygonGate::evaluate(std::vector<bool> &membership)
{
    std::cout << "Applying PolygonGate: " << name;
    enum Orientation { horizontal, vertical, upward, downward };
    struct Edge
    {
        float x0, y0, x1, y1, ymin, ymax, slope;
        Orientation orientation;
        Edge(float X0, float Y0, float X1, float Y1) : x0(X0), y0(Y0), x1(X1), y1(Y1), ymin(std::min(Y0,Y1)), ymax(std::max(Y0,Y1))
        {
            if (x0 > x1)
            {
                std::swap(x0,x1);
                std::swap(y0,y1);
            }
            if (x0 == x1)
                orientation = vertical;
            else if (y0 == y1)
                orientation = horizontal;
            else
            {
                slope = (y1 - y0) / (x1 - x0);
                if (slope > 0)
                    orientation = upward;
                else
                    orientation = downward;
            }
        }
    };
    
    std::vector<Edge> edges;
    if (edges.empty())
    {
        float x0, y0, x1, y1;
        for (size_t i = 0; i < vertices.size() - 1; ++i)
        {
            x0 = dimensions[0].transform->scale(vertices[i][0]);
            y0 = dimensions[1].transform->scale(vertices[i][1]);
            x1 = dimensions[0].transform->scale(vertices[i+1][0]);
            y1 = dimensions[1].transform->scale(vertices[i+1][1]);
            edges.emplace_back(x0, y0, x1, y1);
        }
        x0 = dimensions[0].transform->scale(vertices[vertices.size()-1][0]);
        y0 = dimensions[1].transform->scale(vertices[vertices.size()-1][1]);
        x1 = dimensions[0].transform->scale(vertices[0][0]);
        y1 = dimensions[1].transform->scale(vertices[0][1]);
        edges.emplace_back(x0, y0, x1, y1);

        std::sort(edges.begin(), edges.end(), [](const Edge &a, const Edge & b) -> bool
        { return a.x0 < b.x0; });
    }

    float *xdata = data[0].get()->data();
    float *ydata = data[1].get()->data();
    size_t count = 0;
    for (size_t i = 0; i < membership.size(); ++i)
    {
        if (!membership[i])
            continue;
        float x = xdata[i];
        float y = ydata[i];
        unsigned winding = 0;
        for (const Edge &edge : edges)
        {
            if (edge.x0 > x)
                break;  // nothing left in the list is possible
            if (x > edge.x1)
            {
                if (y >= edge.ymin && y < edge.ymax)
                    ++winding;
            }
            else
            {
                switch (edge.orientation)
                {
                case horizontal:
                    if (y == edge.y0)
                        ++winding;
                    break;
                case vertical:
                    if (y >= edge.ymin && y < edge.ymax)
                        ++winding;
                    break;
                case upward:
                    if (y >= edge.ymin && (y - edge.y0) <= edge.slope * (x - edge.x0))
                        ++winding;
                    break;
                case downward:
                    if (y <= edge.ymax && (y - edge.y0) >= edge.slope * (x - edge.x0))
                        ++winding;
                    break;
                }
            }
        }
        if ((winding & 1) == 0)
            membership[i] = false;
        else
            ++count;
    }
    std::cout << " found " << count << " events." << std::endl;
}

BooleanGate::~BooleanGate() = default;

void BooleanGate::evaluate(std::vector<bool> &membership)
{
    std::cout << "Applying BooleanGate: " << name << std::endl;
    // TODO: Parse and process relevant parameters (logic operations: and, or, not)
}

EllipsoidGate::~EllipsoidGate() = default;

void EllipsoidGate::evaluate(std::vector<bool> &membership)
{
    std::cout << "Applying EllipsoidGate: " << name << std::endl;
    // TODO: Parse and process relevant parameters (dimensions, focus, distance, etc.)
}

QuadrantGate::~QuadrantGate() = default;

void QuadrantGate::evaluate(std::vector<bool> &membership)
{
    std::cout << "Applying QuadrantGate: " << name << std::endl;
    // TODO: Parse and process relevant parameters (dimensions, dividers)
}

// Preorder walk of the gating hierarchy
void walkGatingTree(const std::shared_ptr<Gate> &node, int depth)
{
    if (!node)
        return;

    // Indentation for visualization based on depth
    std::string indent(depth * 2, ' ');
    std::cout << indent << "- ";

    // 1. Visit the node itself
    std::cout << "Gate: " << node->id << "\n";
    for (const auto &dim : node->dimensions)
    {
        std::cout << indent << "  Dim: " << dim.name
                  << " [Scale: " << (dim.scale.empty() ? "None" : dim.scale)
                  << ", Comp: " << (dim.compensation.empty() ? "None" : dim.compensation) << "]\n";
    }

    // 2. Recursively visit the children
    for (const auto &child : node->children)
    {
        walkGatingTree(child, depth + 1);
    }
}

#ifdef BUILD_GATING_TEST
int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <gating-ml.xml>" << std::endl;
        return 1;
    }

    xmlInitParser();

    xmlDocPtr doc = xmlParseFile(argv[1]);
    if (doc == nullptr)
    {
        std::cerr << "XML parsed with errors." << std::endl;
        return 1;
    }

    std::unordered_map<std::string, std::shared_ptr<Gate>> gates;
    std::vector<std::shared_ptr<Gate>> root_gates;

    xmlXPathContextPtr xpathCtx = xmlXPathNewContext(doc);
    if (xpathCtx == nullptr)
    {
        std::cerr << "Error: unable to create new XPath context" << std::endl;
        xmlFreeDoc(doc);
        return 1;
    }

    std::unordered_map<std::string, std::shared_ptr<Transform>> transforms;

    xmlXPathObjectPtr transObj = xmlXPathEvalExpression((const xmlChar *)"//*[local-name()='logicle' or local-name()='hyperlog' or local-name()='lin' or local-name()='log' or local-name()='fasinh' or local-name()='biexp']", xpathCtx);
    if (transObj && transObj->nodesetval)
    {
        for (int i = 0; i < transObj->nodesetval->nodeNr; ++i)
        {
            xmlNodePtr node = transObj->nodesetval->nodeTab[i];

            std::string id;
            xmlChar *idAttr = xmlGetProp(node, (const xmlChar *)"id");
            if (!idAttr)
                idAttr = xmlGetProp(node, (const xmlChar *)"gating:id");
            if (!idAttr && node->parent)
            {
                idAttr = xmlGetProp(node->parent, (const xmlChar *)"id");
                if (!idAttr)
                    idAttr = xmlGetProp(node->parent, (const xmlChar *)"gating:id");
            }
            if (idAttr)
            {
                id = (const char *)idAttr;
                xmlFree(idAttr);
            }

            if (id.empty())
                continue;

            auto get_double = [](xmlNodePtr n, const char *attr_name, double def)
            {
                xmlChar *attr = xmlGetProp(n, (const xmlChar *)attr_name);
                if (attr)
                {
                    double val = def;
                    try
                    {
                        val = std::stod((char *)attr);
                    }
                    catch (...)
                    {
                    }
                    xmlFree(attr);
                    return val;
                }
                return def;
            };

            std::string name = (const char *)node->name;
            std::shared_ptr<Transform> transform;
            try
            {
                if (name.find("logicle") != std::string::npos)
                {
                    transform = std::make_shared<Logicle>(get_double(node, "T", 262144.0), get_double(node, "W", 0.5), get_double(node, "M", 4.5), get_double(node, "A", 0.0));
                }
                else if (name.find("hyperlog") != std::string::npos)
                {
                    transform = std::make_shared<Hyperlog>(get_double(node, "T", 262144.0), get_double(node, "W", 0.5), get_double(node, "M", 4.5), get_double(node, "A", 0.0));
                }
                else if (name.find("lin") != std::string::npos)
                {
                    transform = std::make_shared<Linear>(get_double(node, "T", 262144.0), get_double(node, "A", 0.0));
                }
                else if (name.find("log") != std::string::npos)
                {
                    transform = std::make_shared<Logarithmic>(get_double(node, "T", 262144.0), get_double(node, "M", 4.5));
                }
                else if (name.find("fasinh") != std::string::npos)
                {
                    transform = std::make_shared<Arcsinh>(get_double(node, "T", 262144.0), get_double(node, "M", 4.5), get_double(node, "A", 0.0));
                }
                else if (name.find("biexp") != std::string::npos)
                {
                    transform = std::make_shared<Logicle>(get_double(node, "T", 262144.0), get_double(node, "W", 0.5), get_double(node, "M", 4.5), get_double(node, "A", 0.0));
                }
            }
            catch (...)
            {
            }

            if (transform)
            {
                transforms[id] = transform;
            }
        }
        xmlXPathFreeObject(transObj);
    }

    // Step 1: Extract all gates across the document
    xmlXPathObjectPtr xpathObj = xmlXPathEvalExpression((const xmlChar *)"//*[contains(local-name(), 'Gate')]", xpathCtx);
    if (xpathObj != nullptr && xpathObj->nodesetval != nullptr)
    {
        int size = xpathObj->nodesetval->nodeNr;
        for (int i = 0; i < size; ++i)
        {
            xmlNodePtr node = xpathObj->nodesetval->nodeTab[i];
            std::string name = (const char *)node->name;

            std::string id;
            xmlChar *idAttr = xmlGetProp(node, (const xmlChar *)"gating:id");
            if (!idAttr)
                idAttr = xmlGetProp(node, (const xmlChar *)"id");
            if (!idAttr)
            {
                for (xmlAttrPtr attr = node->properties; attr != nullptr; attr = attr->next)
                {
                    if (xmlStrEqual(attr->name, (const xmlChar *)"id"))
                    {
                        idAttr = xmlNodeListGetString(node->doc, attr->children, 1);
                        break;
                    }
                }
            }
            if (idAttr)
            {
                id = (const char *)idAttr;
                xmlFree(idAttr);
            }

            std::string parent_id;
            xmlChar *parentIdAttr = xmlGetProp(node, (const xmlChar *)"gating:parent_id");
            if (!parentIdAttr)
                parentIdAttr = xmlGetProp(node, (const xmlChar *)"parent_id");
            if (!parentIdAttr)
            {
                for (xmlAttrPtr attr = node->properties; attr != nullptr; attr = attr->next)
                {
                    if (xmlStrEqual(attr->name, (const xmlChar *)"parent_id"))
                    {
                        parentIdAttr = xmlNodeListGetString(node->doc, attr->children, 1);
                        break;
                    }
                }
            }
            if (parentIdAttr)
            {
                parent_id = (const char *)parentIdAttr;
                xmlFree(parentIdAttr);
            }

            std::shared_ptr<Gate> gate;
            if (name.find("RectangleGate") != std::string::npos)
                gate = std::make_shared<RectangleGate>(id, parent_id);
            else if (name.find("PolygonGate") != std::string::npos)
                gate = std::make_shared<PolygonGate>(id, parent_id);
            else if (name.find("BooleanGate") != std::string::npos)
                gate = std::make_shared<BooleanGate>(id, parent_id);
            else if (name.find("EllipsoidGate") != std::string::npos)
                gate = std::make_shared<EllipsoidGate>(id, parent_id);
            else if (name.find("QuadrantGate") != std::string::npos)
                gate = std::make_shared<QuadrantGate>(id, parent_id);

            if (gate && !id.empty())
            {
                xmlXPathObjectPtr dimObj = xmlXPathNodeEval(node, (const xmlChar *)".//*[contains(local-name(), 'fcs-dimension')]", xpathCtx);
                if (dimObj && dimObj->nodesetval)
                {
                    for (int j = 0; j < dimObj->nodesetval->nodeNr; ++j)
                    {
                        xmlNodePtr dimNode = dimObj->nodesetval->nodeTab[j];
                        Gate::Dimension dim;

                        xmlChar *nAttr = xmlGetProp(dimNode, (const xmlChar *)"data-type:name");
                        if (!nAttr)
                            nAttr = xmlGetProp(dimNode, (const xmlChar *)"name");
                        if (nAttr)
                        {
                            dim.name = (char *)nAttr;
                            xmlFree(nAttr);
                        }

                        xmlChar *sAttr = xmlGetProp(dimNode, (const xmlChar *)"data-type:transformation-ref");
                        if (!sAttr)
                            sAttr = xmlGetProp(dimNode, (const xmlChar *)"transformation-ref");
                        if (sAttr)
                        {
                            dim.scale = (char *)sAttr;
                            auto it = transforms.find(dim.scale);
                            if (it != transforms.end())
                            {
                                dim.transform = it->second;
                            }
                            xmlFree(sAttr);
                        }

                        xmlChar *cAttr = xmlGetProp(dimNode, (const xmlChar *)"data-type:compensation-ref");
                        if (!cAttr)
                            cAttr = xmlGetProp(dimNode, (const xmlChar *)"compensation-ref");
                        if (cAttr)
                        {
                            dim.compensation = (char *)cAttr;
                            xmlFree(cAttr);
                        }

                        if (!dim.name.empty())
                        {
                            gate->dimensions.push_back(dim);
                        }
                    }
                }
                if (dimObj)
                    xmlXPathFreeObject(dimObj);

                gates[id] = gate;
            }
        }
        xmlXPathFreeObject(xpathObj);
    }

    xmlXPathFreeContext(xpathCtx);
    xmlFreeDoc(doc);
    xmlCleanupParser();

    // Step 2: Build the hierarchy tree
    for (auto &pair : gates)
    {
        auto &gate = pair.second;
        if (gate->parent_id.empty())
        {
            root_gates.push_back(gate);
        }
        else
        {
            auto parent_it = gates.find(gate->parent_id);
            if (parent_it != gates.end())
                parent_it->second->children.push_back(gate);
            else
            {
                std::cerr << "Warning: Parent ID '" << gate->parent_id << "' not found. Re-rooting gate '" << gate->id << "'\n";
                root_gates.push_back(gate);
            }
        }
    }

    // Step 3: Perform preorder walk starting at the roots
    std::cout << "Starting Preorder Walk of Gating Hierarchy..." << std::endl;
    for (const auto &root : root_gates)
    {
        walkGatingTree(root);
    }

    return 0;
}
#endif