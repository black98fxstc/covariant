#include "FlowJo.hpp"

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

std::vector<std::string> analysis_choices = {"EPP", "Laplace"};

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

            auto it = std::find(sd.populations.begin(), sd.populations.end(), "All");
            if (it != sd.populations.end()) {
                sd.populations.erase(it);
            }
            sd.populations.insert(sd.populations.begin(), "All");

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
    ws.all_populations.erase(std::unique(ws.all_populations.begin(), ws.all_populations.end()), ws.all_populations.end());
    auto it = std::find(ws.all_populations.begin(), ws.all_populations.end(), "All");
    if (it != ws.all_populations.end()) {
        ws.all_populations.erase(it);
    }
    ws.all_populations.insert(ws.all_populations.begin(), "All");

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

void add_laplace_derived_parameter(const std::string &filename, const std::string &sample_name)
{
    xmlDocPtr doc = xmlParseFile(filename.c_str());
    if (!doc)
    {
        std::cerr << "Error: Could not open workspace to add DerivedParameter: " << filename << "\n";
        return;
    }

    xmlXPathContextPtr xpathCtx = xmlXPathNewContext(doc);
    if (!xpathCtx)
    {
        xmlFreeDoc(doc);
        return;
    }

    std::string expr = "//SampleList/Sample[SampleNode/@name='" + sample_name + ".fcs' or SampleNode/@name='" + sample_name + ".FCS' or SampleNode/@name='" + sample_name + "' or .//Keyword[@name='$FIL' and (@value='" + sample_name + ".fcs' or @value='" + sample_name + ".FCS' or @value='" + sample_name + "')]]";
    xmlXPathObjectPtr sampleObj = xmlXPathEvalExpression((const xmlChar *)expr.c_str(), xpathCtx);
    if (sampleObj && sampleObj->nodesetval && sampleObj->nodesetval->nodeNr > 0)
    {
        xmlNodePtr sampleNode = sampleObj->nodesetval->nodeTab[0];

        xmlXPathObjectPtr datasetObj = xmlXPathNodeEval(sampleNode, (const xmlChar *)"./DataSet", xpathCtx);
        std::string uri;
        if (datasetObj && datasetObj->nodesetval && datasetObj->nodesetval->nodeNr > 0)
        {
            xmlChar *uriAttr = xmlGetProp(datasetObj->nodesetval->nodeTab[0], (const xmlChar *)"uri");
            if (uriAttr)
            {
                uri = (char *)uriAttr;
                xmlFree(uriAttr);
            }
        }
        if (datasetObj) xmlXPathFreeObject(datasetObj);

        size_t pos = uri.find(".fcs");
        if (pos == std::string::npos) pos = uri.find(".FCS");
        if (pos != std::string::npos)
        {
            uri.replace(pos, 4, ".csv");
        }
        
        size_t last_slash = uri.find_last_of('/');
        if (last_slash != std::string::npos) {
            std::string filename_part = uri.substr(last_slash + 1);
            size_t p;
            while ((p = filename_part.find("%20")) != std::string::npos) filename_part.replace(p, 3, "_");
            while ((p = filename_part.find(" ")) != std::string::npos) filename_part.replace(p, 1, "_");
            uri = uri.substr(0, last_slash + 1) + filename_part;
        } else {
            size_t p;
            while ((p = uri.find("%20")) != std::string::npos) uri.replace(p, 3, "_");
            while ((p = uri.find(" ")) != std::string::npos) uri.replace(p, 1, "_");
        }

        xmlXPathObjectPtr dpObj = xmlXPathNodeEval(sampleNode, (const xmlChar *)"./DerivedParameters", xpathCtx);
        xmlNodePtr dpNode = nullptr;
        if (dpObj && dpObj->nodesetval && dpObj->nodesetval->nodeNr > 0)
        {
            dpNode = dpObj->nodesetval->nodeTab[0];
        }
        else
        {
            dpNode = xmlNewNode(nullptr, (const xmlChar *)"DerivedParameters");
            xmlNodePtr keywordsNode = nullptr;
            for (xmlNodePtr child = sampleNode->children; child; child = child->next)
            {
                if (child->type == XML_ELEMENT_NODE && (xmlStrcmp(child->name, (const xmlChar *)"Keywords") == 0 || xmlStrcmp(child->name, (const xmlChar *)"SampleNode") == 0))
                {
                    keywordsNode = child;
                    break;
                }
            }
            if (keywordsNode)
            {
                xmlAddPrevSibling(keywordsNode, dpNode);
                xmlAddPrevSibling(keywordsNode, xmlNewText((const xmlChar *)"\n       "));
            }
            else
            {
                xmlAddChild(sampleNode, dpNode);
            }
        }
        if (dpObj) xmlXPathFreeObject(dpObj);

        bool modified = false;

        xmlXPathObjectPtr pObj = xmlXPathNodeEval(dpNode, (const xmlChar *)"./DerivedParameter[@name='Laplace']", xpathCtx);
        if (!pObj || !pObj->nodesetval || pObj->nodesetval->nodeNr == 0)
        {
            xmlNodePtr paramNode = xmlNewNode(nullptr, (const xmlChar *)"DerivedParameter");
            xmlSetProp(paramNode, (const xmlChar *)"name", (const xmlChar *)"Laplace");
            xmlSetProp(paramNode, (const xmlChar *)"type", (const xmlChar *)"importCsv");
            xmlSetProp(paramNode, (const xmlChar *)"importFile", (const xmlChar *)uri.c_str());
            xmlSetProp(paramNode, (const xmlChar *)"range", (const xmlChar *)"1024");
            xmlSetProp(paramNode, (const xmlChar *)"columnIndex", (const xmlChar *)"1");

            xmlNodePtr transNode = xmlNewNode(nullptr, (const xmlChar *)"Transform");

            xmlNsPtr transformsNs = xmlSearchNs(doc, doc->children, (const xmlChar *)"transforms");
            xmlNsPtr datatypeNs = xmlSearchNs(doc, doc->children, (const xmlChar *)"data-type");

            xmlNodePtr linearNode = nullptr;
            if (transformsNs) {
                linearNode = xmlNewNode(transformsNs, (const xmlChar *)"linear");
                xmlSetNsProp(linearNode, transformsNs, (const xmlChar *)"minRange", (const xmlChar *)"0");
                xmlSetNsProp(linearNode, transformsNs, (const xmlChar *)"maxRange", (const xmlChar *)"1024");
            } else {
                linearNode = xmlNewNode(nullptr, (const xmlChar *)"transforms:linear");
                xmlSetProp(linearNode, (const xmlChar *)"transforms:minRange", (const xmlChar *)"0");
                xmlSetProp(linearNode, (const xmlChar *)"transforms:maxRange", (const xmlChar *)"1024");
            }
            xmlSetProp(linearNode, (const xmlChar *)"gain", (const xmlChar *)"1");

            xmlNodePtr pNode = nullptr;
            if (datatypeNs) {
                pNode = xmlNewNode(datatypeNs, (const xmlChar *)"parameter");
                xmlSetNsProp(pNode, datatypeNs, (const xmlChar *)"name", (const xmlChar *)"Laplace");
            } else {
                pNode = xmlNewNode(nullptr, (const xmlChar *)"data-type:parameter");
                xmlSetProp(pNode, (const xmlChar *)"data-type:name", (const xmlChar *)"Laplace");
            }

            xmlAddChild(linearNode, xmlNewText((const xmlChar *)"\n               "));
            xmlAddChild(linearNode, pNode);
            xmlAddChild(linearNode, xmlNewText((const xmlChar *)"\n             "));

            xmlAddChild(transNode, xmlNewText((const xmlChar *)"\n             "));
            xmlAddChild(transNode, linearNode);
            xmlAddChild(transNode, xmlNewText((const xmlChar *)"\n           "));

            xmlAddChild(paramNode, xmlNewText((const xmlChar *)"\n           "));
            xmlAddChild(paramNode, transNode);
            xmlAddChild(paramNode, xmlNewText((const xmlChar *)"\n         "));

            xmlAddChild(dpNode, xmlNewText((const xmlChar *)"\n         "));
            xmlAddChild(dpNode, paramNode);
            xmlAddChild(dpNode, xmlNewText((const xmlChar *)"\n       "));

            modified = true;
        }
        if (pObj) xmlXPathFreeObject(pObj);

        xmlXPathObjectPtr transBlockObj = xmlXPathNodeEval(sampleNode, (const xmlChar *)"./Transformations", xpathCtx);
        xmlNodePtr transBlockNode = nullptr;
        if (transBlockObj && transBlockObj->nodesetval && transBlockObj->nodesetval->nodeNr > 0)
        {
            transBlockNode = transBlockObj->nodesetval->nodeTab[0];
        }
        else
        {
            transBlockNode = xmlNewNode(nullptr, (const xmlChar *)"Transformations");
            xmlNodePtr keywordsNode = nullptr;
            for (xmlNodePtr child = sampleNode->children; child; child = child->next)
            {
                if (child->type == XML_ELEMENT_NODE && (xmlStrcmp(child->name, (const xmlChar *)"Keywords") == 0 || xmlStrcmp(child->name, (const xmlChar *)"DerivedParameters") == 0 || xmlStrcmp(child->name, (const xmlChar *)"SampleNode") == 0))
                {
                    keywordsNode = child;
                    break;
                }
            }
            if (keywordsNode)
            {
                xmlAddPrevSibling(keywordsNode, transBlockNode);
                xmlAddPrevSibling(keywordsNode, xmlNewText((const xmlChar *)"\n       "));
            }
            else
            {
                xmlAddChild(sampleNode, transBlockNode);
            }
        }
        if (transBlockObj) xmlXPathFreeObject(transBlockObj);

        std::string transExpr = "./*[local-name()='linear']/*[local-name()='parameter' and (@name='Laplace' or @data-type:name='Laplace')]";
        xmlXPathObjectPtr tCheckObj = xmlXPathNodeEval(transBlockNode, (const xmlChar *)transExpr.c_str(), xpathCtx);
        if (!tCheckObj || !tCheckObj->nodesetval || tCheckObj->nodesetval->nodeNr == 0)
        {
            xmlNsPtr transformsNs = xmlSearchNs(doc, doc->children, (const xmlChar *)"transforms");
            xmlNsPtr datatypeNs = xmlSearchNs(doc, doc->children, (const xmlChar *)"data-type");

            xmlNodePtr linearNode = nullptr;
            if (transformsNs) {
                linearNode = xmlNewNode(transformsNs, (const xmlChar *)"linear");
                xmlSetNsProp(linearNode, transformsNs, (const xmlChar *)"minRange", (const xmlChar *)"0");
                xmlSetNsProp(linearNode, transformsNs, (const xmlChar *)"maxRange", (const xmlChar *)"1024");
            } else {
                linearNode = xmlNewNode(nullptr, (const xmlChar *)"transforms:linear");
                xmlSetProp(linearNode, (const xmlChar *)"transforms:minRange", (const xmlChar *)"0");
                xmlSetProp(linearNode, (const xmlChar *)"transforms:maxRange", (const xmlChar *)"1024");
            }
            xmlSetProp(linearNode, (const xmlChar *)"gain", (const xmlChar *)"1");

            xmlNodePtr pNode = nullptr;
            if (datatypeNs) {
                pNode = xmlNewNode(datatypeNs, (const xmlChar *)"parameter");
                xmlSetNsProp(pNode, datatypeNs, (const xmlChar *)"name", (const xmlChar *)"Laplace");
            } else {
                pNode = xmlNewNode(nullptr, (const xmlChar *)"data-type:parameter");
                xmlSetProp(pNode, (const xmlChar *)"data-type:name", (const xmlChar *)"Laplace");
            }

            xmlAddChild(linearNode, xmlNewText((const xmlChar *)"\n           "));
            xmlAddChild(linearNode, pNode);
            xmlAddChild(linearNode, xmlNewText((const xmlChar *)"\n         "));

            xmlAddChild(transBlockNode, xmlNewText((const xmlChar *)"\n         "));
            xmlAddChild(transBlockNode, linearNode);
            xmlAddChild(transBlockNode, xmlNewText((const xmlChar *)"\n       "));

            modified = true;
        }
        if (tCheckObj) xmlXPathFreeObject(tCheckObj);

        if (modified)
        {
            xmlSaveFormatFile(filename.c_str(), doc, 1);
        }
    }
    if (sampleObj) xmlXPathFreeObject(sampleObj);

    xmlXPathFreeContext(xpathCtx);
    xmlFreeDoc(doc);
}

void add_laplace_gates(const std::string &filename, const std::string &sample_name, const std::string &parent_pop_name, unsigned clusters_found, size_t laplacian_offset, const std::vector<size_t>& cluster_counts, const std::vector<std::string>& selected_vars)
{
    xmlDocPtr doc = xmlParseFile(filename.c_str());
    if (!doc) {
        std::cerr << "Error: Could not open workspace to add gates: " << filename << "\n";
        return;
    }

    xmlXPathContextPtr xpathCtx = xmlXPathNewContext(doc);
    if (!xpathCtx) {
        xmlFreeDoc(doc);
        return;
    }

    std::string sample_expr = "//SampleList/Sample[SampleNode/@name='" + sample_name + ".fcs' or SampleNode/@name='" + sample_name + ".FCS' or SampleNode/@name='" + sample_name + "' or .//Keyword[@name='$FIL' and (@value='" + sample_name + ".fcs' or @value='" + sample_name + ".FCS' or @value='" + sample_name + "')]]";
    xmlXPathObjectPtr sampleObj = xmlXPathEvalExpression((const xmlChar *)sample_expr.c_str(), xpathCtx);
    if (!sampleObj || !sampleObj->nodesetval || sampleObj->nodesetval->nodeNr == 0) {
        if (sampleObj) xmlXPathFreeObject(sampleObj);
        xmlXPathFreeContext(xpathCtx);
        xmlFreeDoc(doc);
        return;
    }
    
    xmlNodePtr sampleNode = sampleObj->nodesetval->nodeTab[0];

    xmlNodePtr parentPopNode = nullptr;
    if (parent_pop_name == "All") {
        xmlXPathObjectPtr snObj = xmlXPathNodeEval(sampleNode, (const xmlChar *)"./SampleNode", xpathCtx);
        if (snObj && snObj->nodesetval && snObj->nodesetval->nodeNr > 0) {
            parentPopNode = snObj->nodesetval->nodeTab[0];
        }
        if (snObj) xmlXPathFreeObject(snObj);
    } else {
        std::string pop_expr = ".//Population[@name='" + parent_pop_name + "']";
        xmlXPathObjectPtr popObj = xmlXPathNodeEval(sampleNode, (const xmlChar *)pop_expr.c_str(), xpathCtx);
        if (popObj && popObj->nodesetval && popObj->nodesetval->nodeNr > 0) {
            parentPopNode = popObj->nodesetval->nodeTab[0];
        }
        if (popObj) xmlXPathFreeObject(popObj);
    }

    if (!parentPopNode) {
        xmlXPathFreeObject(sampleObj);
        xmlXPathFreeContext(xpathCtx);
        xmlFreeDoc(doc);
        return;
    }

    xmlNodePtr subpopsNode = nullptr;
    xmlNodePtr parentGraphNode = nullptr;
    for (xmlNodePtr child = parentPopNode->children; child; child = child->next) {
        if (child->type == XML_ELEMENT_NODE && xmlStrcmp(child->name, (const xmlChar *)"Subpopulations") == 0) {
            subpopsNode = child;
        } else if (child->type == XML_ELEMENT_NODE && xmlStrcmp(child->name, (const xmlChar *)"Graph") == 0) {
            parentGraphNode = child;
        }
    }
    if (!subpopsNode) {
        subpopsNode = xmlNewNode(nullptr, (const xmlChar *)"Subpopulations");
        xmlAddChild(parentPopNode, xmlNewText((const xmlChar *)"\n           "));
        xmlAddChild(parentPopNode, subpopsNode);
    }

    xmlNsPtr gatingNs = xmlSearchNs(doc, doc->children, (const xmlChar *)"gating");
    xmlNsPtr datatypeNs = xmlSearchNs(doc, doc->children, (const xmlChar *)"data-type");

    for (unsigned k = 0; k <= clusters_found; ++k) {
        size_t klass = k + laplacian_offset;
        
        std::string pop_name = (k == 0) ? "Laplace.Ambiguous" : "Laplace.Cluster" + std::to_string(k);
        
        std::string check_expr = "./Population[@name='" + pop_name + "']";
        xmlXPathObjectPtr pCheckObj = xmlXPathNodeEval(subpopsNode, (const xmlChar*)check_expr.c_str(), xpathCtx);
        bool exists = (pCheckObj && pCheckObj->nodesetval && pCheckObj->nodesetval->nodeNr > 0);
        if (pCheckObj) xmlXPathFreeObject(pCheckObj);
        if (exists) continue;

        xmlNodePtr popNode = xmlNewNode(nullptr, (const xmlChar *)"Population");
        xmlSetProp(popNode, (const xmlChar *)"name", (const xmlChar *)pop_name.c_str());
        xmlSetProp(popNode, (const xmlChar *)"expanded", (const xmlChar *)"1");
        std::string count_str = k < cluster_counts.size() ? std::to_string(cluster_counts[k]) : "0";
        xmlSetProp(popNode, (const xmlChar *)"count", (const xmlChar *)count_str.c_str());

        if (parentGraphNode) {
            xmlAddChild(popNode, xmlNewText((const xmlChar *)"\n             "));
            xmlNodePtr copiedGraphNode = xmlCopyNode(parentGraphNode, 1);
            xmlAddChild(popNode, copiedGraphNode);
            
            if (selected_vars.size() >= 2) {
                for (xmlNodePtr child = copiedGraphNode->children; child; child = child->next) {
                    if (child->type == XML_ELEMENT_NODE && xmlStrcmp(child->name, (const xmlChar *)"Axis") == 0) {
                        xmlChar* dimAttr = xmlGetProp(child, (const xmlChar*)"dimension");
                        if (dimAttr) {
                            if (xmlStrcmp(dimAttr, (const xmlChar*)"x") == 0) {
                                xmlSetProp(child, (const xmlChar*)"name", (const xmlChar*)selected_vars[0].c_str());
                            } else if (xmlStrcmp(dimAttr, (const xmlChar*)"y") == 0) {
                                xmlSetProp(child, (const xmlChar*)"name", (const xmlChar*)selected_vars[1].c_str());
                            }
                            xmlFree(dimAttr);
                        }
                    }
                }
            }
        }

        xmlNodePtr gateNode = xmlNewNode(nullptr, (const xmlChar *)"Gate");
        std::string gateId = "LaplaceGate_" + sample_name + "_" + parent_pop_name + "_" + std::to_string(klass);
        xmlSetProp(gateNode, (const xmlChar *)"gating:id", (const xmlChar *)gateId.c_str());

        xmlNodePtr rectGateNode = xmlNewNode(gatingNs, (const xmlChar *)"RectangleGate");
        xmlSetProp(rectGateNode, (const xmlChar *)"eventsInside", (const xmlChar *)"1");
        
        xmlNodePtr dimNode = xmlNewNode(gatingNs, (const xmlChar *)"dimension");
        std::string min_val = std::to_string(6 + 16 * klass);
        std::string max_val = std::to_string(10 + 16 * klass);
        xmlSetProp(dimNode, (const xmlChar *)"gating:min", (const xmlChar *)min_val.c_str());
        xmlSetProp(dimNode, (const xmlChar *)"gating:max", (const xmlChar *)max_val.c_str());

        xmlNodePtr fcsDimNode = xmlNewNode(datatypeNs, (const xmlChar *)"fcs-dimension");
        xmlSetProp(fcsDimNode, (const xmlChar *)"data-type:name", (const xmlChar *)"Laplace");

        xmlAddChild(dimNode, fcsDimNode);
        xmlAddChild(rectGateNode, dimNode);
        xmlAddChild(gateNode, rectGateNode);
        xmlAddChild(popNode, gateNode);
        
        xmlAddChild(subpopsNode, xmlNewText((const xmlChar *)"\n             "));
        xmlAddChild(subpopsNode, popNode);
    }
    xmlAddChild(subpopsNode, xmlNewText((const xmlChar *)"\n           "));

    xmlSaveFormatFile(filename.c_str(), doc, 1);

    if (sampleObj) xmlXPathFreeObject(sampleObj);
    xmlXPathFreeContext(xpathCtx);
    xmlFreeDoc(doc);
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
        std::string pop_name = ws.all_populations[i];
        pop_container->Add(Maybe(cb, [&, pop_name]
                                 {
            if (selected_samples.empty()) return false;
            return std::find(selected_samples.front()->populations.begin(), selected_samples.front()->populations.end(), pop_name) != selected_samples.front()->populations.end(); }));
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
    auto input_max_clusters = Input(&state.max_clusters_str, "15");
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

        auto sample_win = window(text(" Samples "), sample_container->Render() | vscroll_indicator | frame) | size(WIDTH, EQUAL, 42);
        
        auto combined_content = hbox({
            var_container->Render() | size(WIDTH, EQUAL, 20),
            separator(),
            stain_element | size(WIDTH, EQUAL, 21)
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
            ? window(text(" Populations "), pop_container->Render() | vscroll_indicator | frame) | size(WIDTH, EQUAL, 32)
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
            state.max_clusters = 15;
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