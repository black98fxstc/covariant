#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cstdlib>
#include <filesystem>
#include <algorithm>
#include <map>
#include <cctype>
#include <functional>
#ifndef _WIN32
#include <termios.h>
#include <unistd.h>
#endif

// XML / XSLT headers for future content generation
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <libxslt/xslt.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

#include <ftxui/component/component.hpp>
#include <ftxui/component/screen_interactive.hpp>
#include <ftxui/component/event.hpp>
#include <ftxui/dom/elements.hpp>

#include <Eigen/Dense>

#include "Covariant.hpp"
#include "Samples.hpp"
#include "Gating.hpp"

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
    std::vector<std::string> variables;
    std::vector<std::string> stains;
    std::vector<std::string> populations;
    std::vector<std::shared_ptr<Gate>> gates;
    std::unordered_map<std::string, std::shared_ptr<Transform>> transforms;
    SpilloverMatrix spillover_matrix;
    std::unordered_map<std::string, std::string> gate_id_to_pop;
    std::unordered_map<std::string, std::string> pop_to_gate_id;
    bool selected = false;
    bool enabled = true;
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
    std::vector<SampleData *> samples;
    std::vector<std::string> variables;
    std::vector<std::string> populations;
    int analysis_choice = 0;
    bool cancelled = false;

    std::string smoothing_str = "0.01";
    std::string threshold_str = "0.001";
    std::string max_clusters_str = "50";
    std::string min_events_str = "100";
    float smoothing = 0.01f;
    float threshold = 0.001f;
    unsigned max_clusters = 50;
    size_t min_events = 100;
    unsigned grid_size = 256;
};

const std::vector<std::string> analysis_choices = {
    "Complex sample with many dimensions, analyze with two dimensional methods that can be sorted and visualized.",
    "Sample that cannot be resolved in two dimensions, try higher dimensional methods that cannot be sorted or visualized."};

template <unsigned Dimension>
void do_it(const std::vector<std::vector<float> *> &data, const std::vector<bool> &included, std::vector<unsigned short> &classes, SelectionState &params)
{
    std::vector<size_t> idx;
    for (size_t i = 0; i < included.size(); ++i)
        if (included[i])
            idx.push_back(i);

    Events<Dimension> events;
    events.resize(idx.size());
    for (unsigned d = 0; d < Dimension; ++d)
        for (size_t i = 0; i < events.size(); ++i)
            events[i][d] = (*data[d])[idx[i]];

    std::cout << "Processing " << events.size() << " events..." << std::endl;
    Laplace<Dimension> laplace(params.grid_size);
    for (const auto &e : events)
        laplace.event(e);
    std::cout << "Found " << events.size() << " valid events..." << std::endl;

    std::cout << "Analyzing the sample..." << std::endl;
    laplace.analyze(params.smoothing, params.threshold);
    if (laplace.differentialError() > .0001)
        std::cout << "Differential equation solution is unusually bad " << laplace.differentialError() << std::endl;
    else
        std::cout << "Consistency checks passed..." << std::endl;

    std::cout << "Performing Laplacian clustering..." << std::endl;
    unsigned found = laplace.cluster(params.threshold);
    std::fill(classes.begin(), classes.end(), 0);
    for (size_t i = 0; i < events.size(); ++i)
        classes[idx[i]] = laplace.classify(events[i]);
    std::cout << "Found " << found << " clusters." << std::endl;
};

SpilloverMatrix parse_spillover_matrix(xmlNodePtr matrixNode, xmlXPathContextPtr xpathCtx)
{
    SpilloverMatrix sm;
    xmlChar *idAttr = xmlGetProp(matrixNode, (const xmlChar *)"transforms:id");
    if (!idAttr)
        idAttr = xmlGetProp(matrixNode, (const xmlChar *)"gating:id");
    if (!idAttr)
        idAttr = xmlGetProp(matrixNode, (const xmlChar *)"id");
    if (idAttr)
    {
        sm.id = (char *)idAttr;
        xmlFree(idAttr);
    }

    xmlChar *nameAttr = xmlGetProp(matrixNode, (const xmlChar *)"name");
    if (nameAttr)
    {
        sm.name = (char *)nameAttr;
        xmlFree(nameAttr);
    }

    xmlChar *prefixAttr = xmlGetProp(matrixNode, (const xmlChar *)"prefix");
    if (prefixAttr)
    {
        sm.prefix = (char *)prefixAttr;
        xmlFree(prefixAttr);
    }

    xmlChar *suffixAttr = xmlGetProp(matrixNode, (const xmlChar *)"suffix");
    if (suffixAttr)
    {
        sm.suffix = (char *)suffixAttr;
        xmlFree(suffixAttr);
    }

    xmlXPathObjectPtr paramObj = xmlXPathNodeEval(matrixNode, (const xmlChar *)".//*[local-name()='parameters']/*[local-name()='parameter']", xpathCtx);
    if (paramObj && paramObj->nodesetval)
    {
        for (int i = 0; i < paramObj->nodesetval->nodeNr; ++i)
        {
            xmlNodePtr pNode = paramObj->nodesetval->nodeTab[i];
            std::string pname;
            xmlChar *nAttr = xmlGetProp(pNode, (const xmlChar *)"data-type:name");
            if (!nAttr)
                nAttr = xmlGetProp(pNode, (const xmlChar *)"name");
            if (nAttr)
            {
                pname = (char *)nAttr;
                xmlFree(nAttr);
            }

            std::string pinfix;
            xmlChar *iAttr = xmlGetProp(pNode, (const xmlChar *)"userProvidedCompInfix");
            if (iAttr)
            {
                pinfix = (char *)iAttr;
                xmlFree(iAttr);
            }

            sm.parameters.push_back(pname);
            if (!pname.empty() && !pinfix.empty())
            {
                sm.comp_infix_map[pname] = pinfix;
            }
        }
    }
    if (paramObj)
        xmlXPathFreeObject(paramObj);

    sm.matrix.resize(sm.parameters.size(), std::vector<double>(sm.parameters.size(), 0.0));

    xmlXPathObjectPtr spillObj = xmlXPathNodeEval(matrixNode, (const xmlChar *)".//*[local-name()='spillover']", xpathCtx);
    if (spillObj && spillObj->nodesetval)
    {
        for (int i = 0; i < spillObj->nodesetval->nodeNr; ++i)
        {
            xmlNodePtr sNode = spillObj->nodesetval->nodeTab[i];
            std::string row_param;
            xmlChar *nAttr = xmlGetProp(sNode, (const xmlChar *)"data-type:parameter");
            if (!nAttr)
                nAttr = xmlGetProp(sNode, (const xmlChar *)"parameter");
            if (nAttr)
            {
                row_param = (char *)nAttr;
                xmlFree(nAttr);
            }

            int row_idx = -1;
            for (size_t p = 0; p < sm.parameters.size(); ++p)
            {
                if (sm.parameters[p] == row_param)
                {
                    row_idx = p;
                    break;
                }
            }

            if (row_idx >= 0)
            {
                xmlXPathObjectPtr coefObj = xmlXPathNodeEval(sNode, (const xmlChar *)".//*[local-name()='coefficient']", xpathCtx);
                if (coefObj && coefObj->nodesetval)
                {
                    for (int j = 0; j < coefObj->nodesetval->nodeNr; ++j)
                    {
                        xmlNodePtr cNode = coefObj->nodesetval->nodeTab[j];
                        std::string col_param;
                        xmlChar *cnAttr = xmlGetProp(cNode, (const xmlChar *)"data-type:parameter");
                        if (!cnAttr)
                            cnAttr = xmlGetProp(cNode, (const xmlChar *)"parameter");
                        if (cnAttr)
                        {
                            col_param = (char *)cnAttr;
                            xmlFree(cnAttr);
                        }

                        double val = 0.0;
                        xmlChar *vAttr = xmlGetProp(cNode, (const xmlChar *)"transforms:value");
                        if (!vAttr)
                            vAttr = xmlGetProp(cNode, (const xmlChar *)"value");
                        if (vAttr)
                        {
                            try
                            {
                                val = std::stod((char *)vAttr);
                            }
                            catch (...)
                            {
                            }
                            xmlFree(vAttr);
                        }

                        int col_idx = -1;
                        for (size_t p = 0; p < sm.parameters.size(); ++p)
                        {
                            if (sm.parameters[p] == col_param)
                            {
                                col_idx = p;
                                break;
                            }
                        }
                        if (col_idx >= 0)
                        {
                            sm.matrix[row_idx][col_idx] = val;
                        }
                    }
                }
                if (coefObj)
                    xmlXPathFreeObject(coefObj);
            }
        }
    }
    if (spillObj)
        xmlXPathFreeObject(spillObj);

    return sm;
}

std::string find_workspace(int argc, char *argv[])
{
    std::string filename;
    if (argc >= 2)
    {
        filename = argv[1];
        if (!std::filesystem::exists(filename))
        {
            if (std::filesystem::exists(filename + ".wsp"))
                filename += ".wsp";
        }
    }

    // If not provided or missing, locate the most recent .wsp file
    if (filename.empty() || !std::filesystem::exists(filename))
    {
        std::filesystem::file_time_type latest_time = std::filesystem::file_time_type::min();
        for (const auto &entry : std::filesystem::directory_iterator("."))
        {
            if (entry.is_regular_file() && entry.path().extension() == ".wsp")
            {
                auto mtime = entry.last_write_time();
                if (mtime > latest_time || filename.empty())
                {
                    latest_time = mtime;
                    filename = entry.path().string();
                }
            }
        }
    }

    if (filename.empty() || !std::filesystem::exists(filename))
    {
        std::cerr << "Error: No workspace (.wsp) file found or specified.\n";
        return "";
    }
    return filename;
}

Workspace parse_workspace(const std::string &filename)
{
    Workspace ws;
    ws.filename = filename;

    xmlDocPtr doc = xmlParseFile(filename.c_str());
    if (!doc)
    {
        std::cerr << "Error: Could not parse workspace file: " << filename << "\n";
        return ws;
    }
    xmlXPathContextPtr xpathCtx = xmlXPathNewContext(doc);
    if (!xpathCtx)
    {
        std::cerr << "Error: Could not create XPath context.\n";
        xmlFreeDoc(doc);
        return ws;
    }

    // Fetch Samples or SampleNodes
    xmlXPathObjectPtr samplesObj = xmlXPathEvalExpression((const xmlChar *)"//SampleList/Sample", xpathCtx);
    if (samplesObj && samplesObj->nodesetval)
    {
        for (int i = 0; i < samplesObj->nodesetval->nodeNr; ++i)
        {
            xmlNodePtr sampleNode = samplesObj->nodesetval->nodeTab[i];
            SampleData sd;

            xpathCtx->node = sampleNode;
            xmlXPathObjectPtr transObj = xmlXPathEvalExpression((const xmlChar *)".//*[local-name()='logicle' or local-name()='hyperlog' or local-name()='lin' or local-name()='linear' or local-name()='log' or local-name()='fasinh' or local-name()='biexp' or local-name()='biex']", xpathCtx);
            if (transObj && transObj->nodesetval)
            {
                for (int t = 0; t < transObj->nodesetval->nodeNr; ++t)
                {
                    xmlNodePtr node = transObj->nodesetval->nodeTab[t];

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
                        id = (char *)idAttr;
                        xmlFree(idAttr);
                    }

                    if (id.empty())
                    {
                        // Fallback to checking the parameter name since they don't explicitly carry IDs
                        xmlXPathObjectPtr paramObj = xmlXPathNodeEval(node, (const xmlChar *)".//*[local-name()='parameter']", xpathCtx);
                        if (paramObj && paramObj->nodesetval && paramObj->nodesetval->nodeNr > 0)
                        {
                            xmlNodePtr paramNode = paramObj->nodesetval->nodeTab[0];
                            xmlChar *nameAttr = xmlGetProp(paramNode, (const xmlChar *)"data-type:name");
                            if (!nameAttr)
                                nameAttr = xmlGetProp(paramNode, (const xmlChar *)"name");
                            if (nameAttr)
                            {
                                id = (char *)nameAttr;
                                xmlFree(nameAttr);
                            }
                        }
                        if (paramObj)
                            xmlXPathFreeObject(paramObj);
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
                        else if (name.find("lin") != std::string::npos || name.find("linear") != std::string::npos)
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
                        else if (name.find("biexp") != std::string::npos || name.find("biex") != std::string::npos)
                        {
                            transform = std::make_shared<Logicle>(get_double(node, "T", 262144.0), get_double(node, "W", 0.5), get_double(node, "M", 4.5), get_double(node, "A", 0.0));
                        }
                    }
                    catch (...)
                    {
                    }

                    if (transform)
                    {
                        sd.transforms[id] = transform;
                    }
                }
            }
            if (transObj)
                xmlXPathFreeObject(transObj);

            xmlXPathObjectPtr sampleMatObj = xmlXPathNodeEval(sampleNode, (const xmlChar *)".//*[local-name()='spilloverMatrix']", xpathCtx);
            if (sampleMatObj && sampleMatObj->nodesetval && sampleMatObj->nodesetval->nodeNr > 0)
            {
                sd.spillover_matrix = parse_spillover_matrix(sampleMatObj->nodesetval->nodeTab[0], xpathCtx);
            }
            if (sampleMatObj)
                xmlXPathFreeObject(sampleMatObj);

            // Try fetching name from 'name' attribute
            xmlChar *nameAttr = xmlGetProp(sampleNode, (const xmlChar *)"name");
            if (nameAttr)
            {
                sd.name = (char *)nameAttr;
                xmlFree(nameAttr);
            }
            else
            {
                // Fallback to finding the $FIL keyword
                xpathCtx->node = sampleNode;
                xmlXPathObjectPtr filObj = xmlXPathEvalExpression((const xmlChar *)".//Keyword[@name='$FIL']/@value", xpathCtx);
                if (filObj && filObj->nodesetval && filObj->nodesetval->nodeNr > 0)
                {
                    sd.name = (char *)xmlNodeGetContent(filObj->nodesetval->nodeTab[0]);
                }
                else
                {
                    sd.name = "Sample_" + std::to_string(i);
                }
                if (filObj)
                    xmlXPathFreeObject(filObj);
            }

            size_t fcs_pos = sd.name.find(".fcs");
            if (fcs_pos != std::string::npos)
                sd.name.erase(fcs_pos, 4);
            size_t FCS_pos = sd.name.find(".FCS");
            if (FCS_pos != std::string::npos)
                sd.name.erase(FCS_pos, 4);

            // Collect $PnN keywords (variables) and $PnS keywords (stains)
            xpathCtx->node = sampleNode;
            xmlXPathObjectPtr varObj = xmlXPathEvalExpression((const xmlChar *)".//Keyword[starts-with(@name, '$P') and substring(@name, string-length(@name)) = 'N']", xpathCtx);
            if (varObj && varObj->nodesetval)
            {
                for (int j = 0; j < varObj->nodesetval->nodeNr; ++j)
                {
                    xmlNodePtr kwNode = varObj->nodesetval->nodeTab[j];
                    xmlChar *nameAttr = xmlGetProp(kwNode, (const xmlChar *)"name");
                    xmlChar *valAttr = xmlGetProp(kwNode, (const xmlChar *)"value");
                    if (nameAttr && valAttr)
                    {
                        std::string varName = (char *)nameAttr;
                        sd.variables.push_back((char *)valAttr);

                        std::string num = varName.substr(2, varName.length() - 3);
                        std::string stainPath = ".//Keyword[@name='$P" + num + "S']/@value";
                        xmlXPathObjectPtr sObj = xmlXPathEvalExpression((const xmlChar *)stainPath.c_str(), xpathCtx);
                        if (sObj && sObj->nodesetval && sObj->nodesetval->nodeNr > 0)
                        {
                            xmlChar *sVal = xmlNodeGetContent(sObj->nodesetval->nodeTab[0]);
                            if (sVal)
                            {
                                sd.stains.push_back((char *)sVal);
                                xmlFree(sVal);
                            }
                            else
                            {
                                sd.stains.push_back("");
                            }
                        }
                        else
                        {
                            sd.stains.push_back("");
                        }
                        if (sObj)
                            xmlXPathFreeObject(sObj);
                    }
                    if (nameAttr)
                        xmlFree(nameAttr);
                    if (valAttr)
                        xmlFree(valAttr);
                }
            }
            if (varObj)
                xmlXPathFreeObject(varObj);

            // Collect Populations
            xmlXPathObjectPtr popObj = xmlXPathEvalExpression((const xmlChar *)".//Population", xpathCtx);
            if (popObj && popObj->nodesetval)
            {
                for (int j = 0; j < popObj->nodesetval->nodeNr; ++j)
                {
                    xmlChar *popName = xmlGetProp(popObj->nodesetval->nodeTab[j], (const xmlChar *)"name");
                    if (popName)
                    {
                        std::string pname = (char *)popName;
                        pname.erase(std::remove_if(pname.begin(), pname.end(), [](unsigned char c)
                                                   { return c < 32 || c == 127; }),
                                    pname.end());
                        auto start = pname.find_first_not_of(" ");
                        if (start != std::string::npos)
                        {
                            pname = pname.substr(start);
                        }
                        else
                        {
                            pname = "";
                        }
                        auto end = pname.find_last_not_of(" ");
                        if (end != std::string::npos)
                        {
                            pname = pname.substr(0, end + 1);
                        }

                        if (!pname.empty())
                        {
                            if (std::find(sd.populations.begin(), sd.populations.end(), pname) == sd.populations.end())
                            {
                                sd.populations.push_back(pname);
                                std::shared_ptr<Gate> gate;
                                std::string id;
                                std::string parent_id;
                                xpathCtx->node = popObj->nodesetval->nodeTab[j];
                                xmlXPathObjectPtr gateObj = xmlXPathEvalExpression((const xmlChar *)".//*[contains(local-name(), 'Gate')]", xpathCtx);
                                if (gateObj && gateObj->nodesetval)
                                {
                                    for (int k = 0; k < gateObj->nodesetval->nodeNr; ++k)
                                    {
                                        xmlNodePtr gateNode = gateObj->nodesetval->nodeTab[k];
                                        std::string gname = (const char *)gateNode->name;
                                        // get the graph ids
                                        if (gname == "Gate")
                                        {
                                            xmlChar *idAttr = xmlGetProp(gateNode, (const xmlChar *)"gating:id");
                                            if (!idAttr)
                                                idAttr = xmlGetProp(gateNode, (const xmlChar *)"id");
                                            if (idAttr)
                                            {
                                                id = (const char *)idAttr;
                                                xmlFree(idAttr);
                                            }

                                            xmlChar *parentIdAttr = xmlGetProp(gateNode, (const xmlChar *)"gating:parent_id");
                                            if (!parentIdAttr)
                                                parentIdAttr = xmlGetProp(gateNode, (const xmlChar *)"parent_id");
                                            if (parentIdAttr)
                                            {
                                                parent_id = (const char *)parentIdAttr;
                                                xmlFree(parentIdAttr);
                                            }

                                            sd.gate_id_to_pop[id] = pname;
                                            sd.pop_to_gate_id[pname] = id;
                                            continue;
                                        }
                                        // in need the min and max for each dimension
                                        else if (gname.find("RectangleGate") != std::string::npos)
                                            gate = std::make_shared<RectangleGate>(id, parent_id);
                                        // I need the verticies of the polygon
                                        else if (gname.find("PolygonGate") != std::string::npos)
                                        {
                                            auto poly = std::make_shared<PolygonGate>(id, parent_id);
                                            xmlXPathObjectPtr vObj = xmlXPathNodeEval(gateNode, (const xmlChar *)".//*[local-name()='vertex']", xpathCtx);
                                            if (vObj && vObj->nodesetval)
                                            {
                                                for (int v = 0; v < vObj->nodesetval->nodeNr; ++v)
                                                {
                                                    std::vector<double> vertex;
                                                    xmlXPathObjectPtr coordObj = xmlXPathNodeEval(vObj->nodesetval->nodeTab[v], (const xmlChar *)".//*[local-name()='coordinate']", xpathCtx);
                                                    if (coordObj && coordObj->nodesetval)
                                                    {
                                                        for (int c = 0; c < coordObj->nodesetval->nodeNr; ++c)
                                                        {
                                                            xmlChar *valAttr = xmlGetProp(coordObj->nodesetval->nodeTab[c], (const xmlChar *)"data-type:value");
                                                            if (!valAttr)
                                                                valAttr = xmlGetProp(coordObj->nodesetval->nodeTab[c], (const xmlChar *)"value");
                                                            if (valAttr)
                                                            {
                                                                try
                                                                {
                                                                    vertex.push_back(std::stod((char *)valAttr));
                                                                }
                                                                catch (...)
                                                                {
                                                                }
                                                                xmlFree(valAttr);
                                                            }
                                                        }
                                                    }
                                                    if (coordObj)
                                                        xmlXPathFreeObject(coordObj);
                                                    poly->vertices.push_back(vertex);
                                                }
                                            }
                                            if (vObj)
                                                xmlXPathFreeObject(vObj);
                                            gate = poly;
                                        }
                                        // do not need for now
                                        else if (gname.find("BooleanGate") != std::string::npos)
                                            gate = std::make_shared<BooleanGate>(id, parent_id);
                                        // I need the mean and covariance matrix
                                        else if (gname.find("EllipsoidGate") != std::string::npos)
                                        {
                                            auto ellip = std::make_shared<EllipsoidGate>(id, parent_id);
                                            xmlXPathObjectPtr meanObj = xmlXPathNodeEval(gateNode, (const xmlChar *)".//*[local-name()='mean']//*[local-name()='coordinate']", xpathCtx);
                                            if (meanObj && meanObj->nodesetval)
                                            {
                                                for (int c = 0; c < meanObj->nodesetval->nodeNr; ++c)
                                                {
                                                    xmlChar *valAttr = xmlGetProp(meanObj->nodesetval->nodeTab[c], (const xmlChar *)"data-type:value");
                                                    if (!valAttr)
                                                        valAttr = xmlGetProp(meanObj->nodesetval->nodeTab[c], (const xmlChar *)"value");
                                                    if (valAttr)
                                                    {
                                                        try
                                                        {
                                                            ellip->mean.push_back(std::stod((char *)valAttr));
                                                        }
                                                        catch (...)
                                                        {
                                                        }
                                                        xmlFree(valAttr);
                                                    }
                                                }
                                            }
                                            if (meanObj)
                                                xmlXPathFreeObject(meanObj);
                                            xmlXPathObjectPtr rowObj = xmlXPathNodeEval(gateNode, (const xmlChar *)".//*[local-name()='covarianceMatrix']//*[local-name()='row']", xpathCtx);
                                            if (rowObj && rowObj->nodesetval)
                                            {
                                                for (int r = 0; r < rowObj->nodesetval->nodeNr; ++r)
                                                {
                                                    std::vector<double> row;
                                                    xmlXPathObjectPtr entryObj = xmlXPathNodeEval(rowObj->nodesetval->nodeTab[r], (const xmlChar *)".//*[local-name()='entry']", xpathCtx);
                                                    if (entryObj && entryObj->nodesetval)
                                                    {
                                                        for (int e = 0; e < entryObj->nodesetval->nodeNr; ++e)
                                                        {
                                                            xmlChar *valAttr = xmlGetProp(entryObj->nodesetval->nodeTab[e], (const xmlChar *)"data-type:value");
                                                            if (!valAttr)
                                                                valAttr = xmlGetProp(entryObj->nodesetval->nodeTab[e], (const xmlChar *)"value");
                                                            if (valAttr)
                                                            {
                                                                try
                                                                {
                                                                    row.push_back(std::stod((char *)valAttr));
                                                                }
                                                                catch (...)
                                                                {
                                                                }
                                                                xmlFree(valAttr);
                                                            }
                                                        }
                                                    }
                                                    if (entryObj)
                                                        xmlXPathFreeObject(entryObj);
                                                    ellip->covariance_matrix.push_back(row);
                                                }
                                            }
                                            if (rowObj)
                                                xmlXPathFreeObject(rowObj);
                                            gate = ellip;
                                        }
                                        // do not need for now
                                        else if (gname.find("QuadrantGate") != std::string::npos)
                                            gate = std::make_shared<QuadrantGate>(id, parent_id);

                                        if (gate)
                                        {
                                            xmlXPathObjectPtr dimObj = xmlXPathNodeEval(gateNode, (const xmlChar *)".//*[local-name()='dimension' or local-name()='Dimension' or local-name()='fcs-dimension']", xpathCtx);
                                            if (dimObj && dimObj->nodesetval)
                                            {
                                                for (int d = 0; d < dimObj->nodesetval->nodeNr; ++d)
                                                {
                                                    xmlNodePtr dimNode = dimObj->nodesetval->nodeTab[d];
                                                    std::string nodeName = (const char *)dimNode->name;
                                                    if (nodeName.find("fcs-dimension") != std::string::npos && dimNode->parent)
                                                    {
                                                        std::string parentName = (const char *)dimNode->parent->name;
                                                        if (parentName.find("dimension") != std::string::npos || parentName.find("Dimension") != std::string::npos)
                                                        {
                                                            continue;
                                                        }
                                                    }

                                                    Gate::Dimension dim;

                                                    xmlChar *minAttr = xmlGetProp(dimNode, (const xmlChar *)"gating:min");
                                                    if (!minAttr)
                                                        minAttr = xmlGetProp(dimNode, (const xmlChar *)"min");
                                                    if (minAttr)
                                                    {
                                                        try
                                                        {
                                                            dim.min_val = std::stod((char *)minAttr);
                                                        }
                                                        catch (...)
                                                        {
                                                        }
                                                        xmlFree(minAttr);
                                                    }

                                                    xmlChar *maxAttr = xmlGetProp(dimNode, (const xmlChar *)"gating:max");
                                                    if (!maxAttr)
                                                        maxAttr = xmlGetProp(dimNode, (const xmlChar *)"max");
                                                    if (maxAttr)
                                                    {
                                                        try
                                                        {
                                                            dim.max_val = std::stod((char *)maxAttr);
                                                        }
                                                        catch (...)
                                                        {
                                                        }
                                                        xmlFree(maxAttr);
                                                    }

                                                    xmlNodePtr paramNode = dimNode;
                                                    xmlXPathObjectPtr fcsDimObj = xmlXPathNodeEval(dimNode, (const xmlChar *)".//*[local-name()='fcs-dimension']", xpathCtx);
                                                    if (fcsDimObj && fcsDimObj->nodesetval && fcsDimObj->nodesetval->nodeNr > 0)
                                                    {
                                                        paramNode = fcsDimObj->nodesetval->nodeTab[0];
                                                    }
                                                    if (fcsDimObj)
                                                        xmlXPathFreeObject(fcsDimObj);

                                                    xmlChar *nAttr = xmlGetProp(paramNode, (const xmlChar *)"data-type:name");
                                                    if (!nAttr)
                                                        nAttr = xmlGetProp(paramNode, (const xmlChar *)"name");
                                                    if (nAttr)
                                                    {
                                                        dim.name = (char *)nAttr;
                                                        xmlFree(nAttr);
                                                    }

                                                    xmlChar *sAttr = xmlGetProp(paramNode, (const xmlChar *)"data-type:transformation-ref");
                                                    if (!sAttr)
                                                        sAttr = xmlGetProp(paramNode, (const xmlChar *)"transformation-ref");
                                                    if (sAttr)
                                                    {
                                                        dim.scale = (char *)sAttr;
                                                        auto it = sd.transforms.find(dim.scale);
                                                        if (it != sd.transforms.end())
                                                        {
                                                            dim.transform = it->second;
                                                        }
                                                        xmlFree(sAttr);
                                                    }

                                                    xmlChar *cAttr = xmlGetProp(paramNode, (const xmlChar *)"data-type:compensation-ref");
                                                    if (!cAttr)
                                                        cAttr = xmlGetProp(paramNode, (const xmlChar *)"compensation-ref");
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
                                            break;
                                        }
                                    }
                                }
                                if (gateObj)
                                    xmlXPathFreeObject(gateObj);
                                gate->name = pname;
                                sd.gates.push_back(gate);
                            }
                            if (std::find(ws.all_populations.begin(), ws.all_populations.end(), pname) == ws.all_populations.end())
                            {
                                ws.all_populations.push_back(pname);
                            }
                        }
                        xmlFree(popName);
                    }
                }
            }
            if (popObj)
                xmlXPathFreeObject(popObj);

            std::unordered_map<std::string, std::shared_ptr<Gate>> id_to_gate;
            for (auto &gate : sd.gates)
                if (gate)
                    id_to_gate[gate->id] = gate;

            for (auto &gate : sd.gates)
            {
                if (gate && !gate->parent_id.empty())
                {
                    auto parent_it = id_to_gate.find(gate->parent_id);
                    if (parent_it != id_to_gate.end())
                    {
                        parent_it->second->children.push_back(gate);
                    }
                }
            }

            ws.samples.push_back(sd);
        }
    }
    if (samplesObj)
        xmlXPathFreeObject(samplesObj);
    xmlXPathFreeContext(xpathCtx);
    xmlFreeDoc(doc);
    xmlCleanupParser();

    std::sort(ws.all_populations.begin(), ws.all_populations.end());

    // Consolidate unique variables across all samples
    for (const auto &s : ws.samples)
    {
        for (const auto &v : s.variables)
        {
            if (std::find(ws.all_variables.begin(), ws.all_variables.end(), v) == ws.all_variables.end())
            {
                ws.all_variables.push_back(v);
            }
        }
    }
    std::sort(ws.all_variables.begin(), ws.all_variables.end());

    return ws;
}

SelectionState build_ftxui_interface(Workspace &ws)
{
    using namespace ftxui;
    SelectionState state;

    auto var_states = std::make_unique<bool[]>(ws.all_variables.size());
    for (size_t i = 0; i < ws.all_variables.size(); ++i)
        var_states[i] = false;

    auto pop_states = std::make_unique<bool[]>(ws.all_populations.size());
    for (size_t i = 0; i < ws.all_populations.size(); ++i)
        pop_states[i] = false;

    int num_samples_selected = 0;
    std::vector<SampleData *> selected_samples;
    int num_vars_selected = 0;
    int num_pops_selected = 0;
    std::vector<std::string> selected_pops;

    auto sample_container = Container::Vertical({});
    for (size_t i = 0; i < ws.samples.size(); ++i)
    {
        auto cb = Checkbox(&ws.samples[i].name, &ws.samples[i].selected);
        auto cb_styled = Renderer(cb, [cb, &ws, i]
                                  {
            auto el = cb->Render();
            if (!ws.samples[i].enabled) el = el | color(Color::GrayDark);
            return el; });
        auto cb_handled = CatchEvent(cb_styled, [&ws, i](ftxui::Event e)
                                     {
            if (!ws.samples[i].enabled) {
                if (e == ftxui::Event::Character(' ') || e == ftxui::Event::Return) return true; // block interaction
            }
            return false; });
        sample_container->Add(cb_handled);
    }

    auto var_container = Container::Vertical({});
    for (size_t i = 0; i < ws.all_variables.size(); ++i)
    {
        auto cb = Checkbox(&ws.all_variables[i], &var_states[i]);
        var_container->Add(Maybe(cb, [&, i]
                                 {
            if (selected_samples.empty()) return false;
            return std::find(selected_samples.front()->variables.begin(), selected_samples.front()->variables.end(), ws.all_variables[i]) != selected_samples.front()->variables.end(); }));
    }

    auto pop_container = Container::Vertical({});
    for (size_t i = 0; i < ws.all_populations.size(); ++i)
    {
        auto cb = Checkbox(&ws.all_populations[i], &pop_states[i]);
        pop_container->Add(Maybe(cb, [&, i]
                                 {
            if (selected_samples.empty()) return false;
            return std::find(selected_samples.front()->populations.begin(), selected_samples.front()->populations.end(), ws.all_populations[i]) != selected_samples.front()->populations.end(); }));
    }

    int analysis_choice = 0;
    // use global analysis_choices vector
    auto choice_container = Radiobox(&analysis_choices, &analysis_choice);
    auto choice_handled = CatchEvent(choice_container, [&](ftxui::Event e)
                                     {
        if (num_vars_selected > 4 && analysis_choice == 0) {
            if (e == ftxui::Event::ArrowDown || e == ftxui::Event::Character('j')) return true; // block switching to 2nd option
        }
        return false; });

    auto input_smoothing = Input(&state.smoothing_str, "0.01");
    auto input_threshold = Input(&state.threshold_str, "0.001");
    auto input_max_clusters = Input(&state.max_clusters_str, "50");
    auto input_min_events = Input(&state.min_events_str, "100");

    auto settings_container = Container::Vertical({input_smoothing,
                                                   input_threshold,
                                                   input_max_clusters,
                                                   input_min_events});

    auto main_layout = Container::Horizontal({sample_container,
                                              Maybe(var_container, [&]
                                                    { return num_samples_selected > 0; }),
                                              Maybe(pop_container, [&]
                                                    { return num_vars_selected >= 2; })});

    auto bottom_container = Container::Horizontal({choice_handled,
                                                   settings_container});

    auto top_level = Container::Vertical({main_layout,
                                          bottom_container});

    auto top_level_handled = CatchEvent(top_level, [&, main_layout, bottom_container](ftxui::Event e)
                                        {
        if (e == ftxui::Event::Tab) {
            if (main_layout->Focused()) {
                bottom_container->TakeFocus();
            } else {
                main_layout->TakeFocus();
            }
            return true;
        }
        return false; });

    auto renderer = Renderer(top_level_handled, [&]
                             {
        // Pre-render state resolution 
        num_samples_selected = 0;
        selected_samples.clear();
        for (auto& s : ws.samples) {
            if (s.selected) {
                num_samples_selected++;
                selected_samples.push_back(&s);
            }
        }
        
        num_vars_selected = 0;
        for (size_t i = 0; i < ws.all_variables.size(); ++i) {
            if (var_states[i] && !selected_samples.empty() && 
                std::find(selected_samples.front()->variables.begin(), selected_samples.front()->variables.end(), ws.all_variables[i]) != selected_samples.front()->variables.end()) {
                num_vars_selected++;
            } else {
                var_states[i] = false; // ensure variables hidden by a new sample are deselected
            }
        }
        
        num_pops_selected = 0;
        selected_pops.clear();
        for (size_t i = 0; i < ws.all_populations.size(); ++i) {
            bool in_sample = !selected_samples.empty() && std::find(selected_samples.front()->populations.begin(), selected_samples.front()->populations.end(), ws.all_populations[i]) != selected_samples.front()->populations.end();
            if (pop_states[i] && in_sample) {
                num_pops_selected++;
                selected_pops.push_back(ws.all_populations[i]);
            } else {
                pop_states[i] = false;
            }
        }
        
        if (num_vars_selected > 4) analysis_choice = 0;
        
        for (auto& s : ws.samples) {
            if (num_pops_selected > 0) {
                bool has_all = true;
                for (const auto& p : selected_pops) {
                    if (std::find(s.populations.begin(), s.populations.end(), p) == s.populations.end()) {
                        has_all = false; break;
                    }
                }
                s.enabled = has_all;
                if (!s.enabled) s.selected = false;
            } else {
                s.enabled = true;
            }
        }

        Elements stain_elements;
        if (!selected_samples.empty()) {
            for (size_t i = 0; i < ws.all_variables.size(); ++i) {
                auto it = std::find(selected_samples.front()->variables.begin(), selected_samples.front()->variables.end(), ws.all_variables[i]);
                if (it != selected_samples.front()->variables.end()) {
                    int idx = std::distance(selected_samples.front()->variables.begin(), it);
                    stain_elements.push_back(text(selected_samples.front()->stains[idx]) | size(HEIGHT, EQUAL, 1));
                }
            }
        }
        auto stain_element = vbox(std::move(stain_elements));

        auto sample_win = window(text(" Samples "), sample_container->Render() | vscroll_indicator | frame) | size(WIDTH, EQUAL, 35);
        
        auto combined_content = hbox({
            var_container->Render() | size(WIDTH, EQUAL, 20),
            separator(),
            stain_element | size(WIDTH, EQUAL, 20)
        });

        auto var_win = num_samples_selected > 0 
            ? window(text(" Detectors            Stains "), combined_content | vscroll_indicator | frame)
            : emptyElement();

        auto settings_win = window(text(" Settings "), vbox({
            hbox({text("Smoothing:    "), input_smoothing->Render() | size(WIDTH, EQUAL, 10)}),
            hbox({text("Threshold:    "), input_threshold->Render() | size(WIDTH, EQUAL, 10)}),
            hbox({text("Max Clusters: "), input_max_clusters->Render() | size(WIDTH, EQUAL, 10)}),
            hbox({text("Min Events:   "), input_min_events->Render() | size(WIDTH, EQUAL, 10)})
        }));
            
        auto pop_win = num_vars_selected >= 2 
            ? window(text(" Populations "), pop_container->Render() | vscroll_indicator | frame) | size(WIDTH, EQUAL, 30)
            : emptyElement();
            
        return vbox({
            text(" Leonard - " + std::filesystem::path(ws.filename).stem().string() + " ") | bold | hcenter,
            separator(),
            hbox({ sample_win, var_win, pop_win }) | size(HEIGHT, LESS_THAN, 20),
            separator(),
            hbox({
                window(text(" Analysis Method "), choice_handled->Render()) | flex,
                settings_win
            }),
            separator(),
            text(" Space: Select | Arrows: Navigate | Tab: Switch Section | Enter: Confirm | Esc: Cancel ") | hcenter
        }) | border; });

    {
        auto screen = ScreenInteractive::TerminalOutput();
        auto main_container = CatchEvent(renderer, [&](ftxui::Event event)
                                         {
            if (event == ftxui::Event::Return) {
                screen.ExitLoopClosure()();
                return true;
            }
            if (event == ftxui::Event::Escape) {
                state.cancelled = true;
                screen.ExitLoopClosure()();
                return true;
            }
            return false; });

        screen.Loop(main_container);
    }

    std::cout << "\x1B[2J\x1B[H";
    std::cout.flush();

#ifndef _WIN32
    // Flush any stray terminal responses (e.g., cursor position) received during the Escape timeout
    tcflush(STDIN_FILENO, TCIFLUSH);
#endif

    if (!state.cancelled)
    {
        for (auto &s : ws.samples)
            if (s.selected)
                state.samples.push_back(&s);
        for (size_t i = 0; i < ws.all_variables.size(); ++i)
        {
            if (var_states[i])
                state.variables.push_back(ws.all_variables[i]);
        }
        for (size_t i = 0; i < ws.all_populations.size(); ++i)
        {
            if (pop_states[i])
                state.populations.push_back(ws.all_populations[i]);
        }
        state.analysis_choice = analysis_choice;

        try
        {
            state.smoothing = std::stof(state.smoothing_str);
        }
        catch (...)
        {
            state.smoothing = 0.01f;
        }
        try
        {
            state.threshold = std::stof(state.threshold_str);
        }
        catch (...)
        {
            state.threshold = 0.001f;
        }
        try
        {
            state.max_clusters = std::stoul(state.max_clusters_str);
        }
        catch (...)
        {
            state.max_clusters = 50;
        }
        try
        {
            state.min_events = std::stoull(state.min_events_str);
        }
        catch (...)
        {
            state.min_events = 100;
        }
    }

    return state;
}

int main(int argc, char *argv[])
{
    std::string filename = find_workspace(argc, argv);
    if (filename.empty())
        return 1;

    Workspace ws = parse_workspace(filename);
    if (ws.samples.empty())
    {
        std::cerr << "Error: No samples found or failed to parse workspace.\n";
        return 1;
    }

    SelectionState selections = build_ftxui_interface(ws);
    if (selections.cancelled)
        return 1;
    DataSet dataset;
    for (const auto *s_ptr : selections.samples)
    {
        const auto &s = *s_ptr;
        std::unordered_map<std::string, std::shared_ptr<Gate>> gate_id_to_gate;
        for (const auto &g : s.gates)
            if (g)
                gate_id_to_gate[g->id] = g;

        dataset.read(s.name);

        // setup spectral unmixing
        std::vector<std::vector<float> *> detector_data;
        for (const auto &v : s.spillover_matrix.parameters)
            detector_data.push_back(dataset.variable[v].data.get());

        const std::vector<std::vector<double>> &spill = s.spillover_matrix.matrix;
        Eigen::MatrixXf spectrum(spill.size(), spill[0].size());
        for (size_t i = 0; i < spill.size(); i++)
            for (size_t j = 0; j < spill[i].size(); j++)
                spectrum(i, j) = spill[i][j];
        Eigen::BDCSVD<Eigen::MatrixXf, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(spectrum);
        Eigen::MatrixXf pseudo(svd.matrixV() * svd.singularValues().asDiagonal().inverse() * svd.matrixU().transpose());
        for (size_t i = 0; i < pseudo.rows(); i++)
        {
            const auto &p = s.spillover_matrix.parameters[i];
            auto v = s.spillover_matrix.comp_infix_map.find(p);
            if (v == s.spillover_matrix.comp_infix_map.end())
                continue;
            Variable &var = dataset.variable[v->second];
            for (size_t j = 0; j < pseudo.cols(); j++)
                var.unmixing.push_back(pseudo(i, j));
        }

        // set the scales on the variables for this sample
        for (auto &var : dataset.variable)
        {
            auto xfrm = s.transforms.find(var.first);
            if (xfrm != s.transforms.end())
                var.second.transform = xfrm->second;
        }

        // cache of transformed data
        std::unordered_map<std::string, std::shared_ptr<std::vector<float>>> transformed_data;

        auto get_transformed = [&](const std::string &var_name) -> std::shared_ptr<std::vector<float>>
        {
            auto it = transformed_data.find(var_name);
            if (it != transformed_data.end())
                return it->second;

            auto &var = dataset.variable[var_name];

            // if the data don't exist we need to do the unmixing first
            if (!var.data && !var.unmixing.empty())
            {
                var.data = std::make_shared<std::vector<float>>(dataset.size(), 0.0f);
                std::vector<float> *bare_data = var.data.get();
                for (size_t i = 0; i < dataset.size(); ++i)
                    for (size_t j = 0; j < detector_data.size(); ++j)
                        (*bare_data)[i] += (*detector_data[j])[i] * var.unmixing[j];
            }

            std::shared_ptr<std::vector<float>> transformed;
            const std::shared_ptr<std::vector<float>> &orig_data = var.data;
            transformed = std::make_shared<std::vector<float>>(orig_data->size());
            for (size_t i = 0; i < orig_data->size(); ++i)
                (*transformed)[i] = var.transform_value((*orig_data)[i]);

            transformed_data[var_name] = transformed; // Cache it for future lookups
            return transformed;
        };

        // generate the appropriate scaled data for analysis
        for (size_t v = 0; v < selections.variables.size(); ++v)
        {
            const auto &vname = selections.variables[v];
            Variable &var = dataset.variable[vname];
            auto transformed = get_transformed(vname);
            if (transformed)
                var.data = transformed;
        }

        // find all the gates and any data they need
        std::vector<std::shared_ptr<Gate>> gates;
        for (size_t i = 0; i < selections.populations.size(); ++i)
        {
            std::string pop_name = selections.populations[i];
            std::string gate_id;
            std::shared_ptr<Gate> gate;
            auto it = s.pop_to_gate_id.find(pop_name);
            if (it != s.pop_to_gate_id.end())
            {
                gate_id = it->second;
                do
                {
                    gate = gate_id_to_gate[gate_id];
                    for (const auto &g : gates)
                        if (g == gate)
                            continue;
                    gates.push_back(gate);
                    for (const auto &dim : gate->dimensions)
                        gate->data.push_back(get_transformed(dim.name));
                    gate_id = gate->parent_id;
                } while (!gate_id.empty());
            }
        }

        // evaluate all the gates
        std::vector<bool> whole(dataset.size(), true);
        for (auto it = gates.rbegin(); it != gates.rend(); ++it)
        {
            auto &gate = *it;
            if (gate->parent_id.empty())
                gate->membership = std::make_shared<std::vector<bool>>(dataset.size(), true);
            else
            {
                std::shared_ptr<Gate> parent = gate_id_to_gate[gate->parent_id];
                gate->membership = std::make_shared<std::vector<bool>>(*(parent->membership.get()));
            }
            gate->evaluate(*gate->membership);
        }

        // someplace for the data to go
        std::vector<unsigned short> *classification = dataset.classification["Laplace"].classifications.get();

        // scaled data marshalled as events for analysis
        int num_vars_selected = selections.variables.size();
        std::vector<std::vector<float> *> data(num_vars_selected);
        for (size_t i = 0; i < num_vars_selected; ++i)
            data[i] = dataset.variable[selections.variables[i]].data.get();

        for (std::string &pop_name : selections.populations)
        {
            std::vector<bool> *subpopulation = dataset.subset[pop_name].membership.get();
            switch (selections.analysis_choice)
            {
            case 0:
                // EPP
                break;
            case 1:
                // Laplace
                switch (num_vars_selected)
                {
                case 2:
                    do_it<2>(data, *subpopulation, *classification, selections);
                    break;
                case 3:
                    do_it<3>(data, *subpopulation, *classification, selections);
                    break;
                case 4:
                    do_it<4>(data, *subpopulation, *classification, selections);
                    break;
                default:
                    // error
                    break;
                }
                break;
            default:
                // error
                break;
            }
        }
    };

    std::cout << "\n=== LEONARD SELECTIONS ===\n";
    if (!selections.samples.empty())
    {
        std::cout << "Samples:";
        for (const auto *s : selections.samples)
            std::cout << " " << s->name;
        std::cout << "\n";
    }
    std::cout << "Variables: ";
    for (const auto &v : selections.variables)
        std::cout << v << " ";

    std::cout << "\nPopulations: ";
    for (const auto &p : selections.populations)
        std::cout << p << " ";

    std::cout << "\n\nSettings:\n";
    std::cout << "  Smoothing: " << selections.smoothing << "\n";
    std::cout << "  Threshold: " << selections.threshold << "\n";
    std::cout << "  Max Clusters: " << selections.max_clusters << "\n";
    std::cout << "  Min Events: " << selections.min_events << "\n";

    std::cout << "\nAnalysis Method:\n"
              << analysis_choices[selections.analysis_choice] << "\n\n";
    return 0;
}
