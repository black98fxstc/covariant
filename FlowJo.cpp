#include "FlowJo.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cstdlib>
#include <filesystem>
#include <algorithm>
#include <map>
#include <cmath>
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
#include <nlohmann/json.hpp>
#include <fstream>

#include <ftxui/component/component.hpp>
#include <ftxui/component/screen_interactive.hpp>
#include <ftxui/component/event.hpp>
#include <ftxui/dom/elements.hpp>

using json = nlohmann::json;

static std::string sanitize_workspace_string(const std::string& value)
{
    std::string sanitized;
    sanitized.reserve(value.size());
    for (unsigned char c : value)
    {
        if (c >= 0x20 && c != 0x7f)
            sanitized.push_back(static_cast<char>(c));
    }
    return sanitized;
}

static std::string xpath_literal(const std::string& value)
{
    if (value.find('\'') == std::string::npos)
        return "'" + value + "'";
    if (value.find('"') == std::string::npos)
        return "\"" + value + "\"";

    std::string expression = "concat(";
    size_t start = 0;
    bool first = true;
    while (start < value.size())
    {
        size_t quote = value.find_first_of("'\"", start);
        std::string part = value.substr(start, quote == std::string::npos ? std::string::npos : quote - start);
        if (!first)
            expression += ",";
        expression += "'" + part + "'";
        first = false;
        if (quote == std::string::npos)
            break;
        expression += value[quote] == '\'' ? "\"'\"" : "'\"'";
        start = quote + 1;
    }
    expression += ")";
    return expression;
}

std::vector<std::string> analysis_choices = {"Exhaustive Projection Pursuit", "Laplacian Clustering", "Covariant Statistics"};

std::string get_settings_path() {
    std::string path;
#ifdef _WIN32
    if (const char* appdata = std::getenv("APPDATA")) {
        path = std::string(appdata) + "\\Leonard";
    }
#else
    if (const char* home = std::getenv("HOME")) {
        path = std::string(home) + "/.config/Leonard";
    }
#endif
    if (!path.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(path, ec);
        path += "/settings.json";
    }
    return path;
}

SpilloverMatrix parse_spillover_matrix(xmlNodePtr matrixNode, xmlXPathContextPtr xpathCtx)
{
    SpilloverMatrix sm;
    xmlChar *idAttr = xmlGetProp(matrixNode, reinterpret_cast<const xmlChar *>("transforms:id"));
    if (!idAttr)
        idAttr = xmlGetProp(matrixNode, reinterpret_cast<const xmlChar *>("gating:id"));
    if (!idAttr)
        idAttr = xmlGetProp(matrixNode, reinterpret_cast<const xmlChar *>("id"));
    if (idAttr)
    {
        sm.id = sanitize_workspace_string(reinterpret_cast<char *>(idAttr));
        xmlFree(idAttr);
    }

    xmlChar *nameAttr = xmlGetProp(matrixNode, reinterpret_cast<const xmlChar *>("name"));
    if (nameAttr)
    {
        sm.name = sanitize_workspace_string(reinterpret_cast<char *>(nameAttr));
        xmlFree(nameAttr);
    }

    xmlChar *prefixAttr = xmlGetProp(matrixNode, reinterpret_cast<const xmlChar *>("prefix"));
    if (prefixAttr)
    {
        sm.prefix = sanitize_workspace_string(reinterpret_cast<char *>(prefixAttr));
        xmlFree(prefixAttr);
    }

    xmlChar *suffixAttr = xmlGetProp(matrixNode, reinterpret_cast<const xmlChar *>("suffix"));
    if (suffixAttr)
    {
        sm.suffix = sanitize_workspace_string(reinterpret_cast<char *>(suffixAttr));
        xmlFree(suffixAttr);
    }

    xmlXPathObjectPtr paramObj = xmlXPathNodeEval(matrixNode, reinterpret_cast<const xmlChar *>(".//*[local-name()='parameters']/*[local-name()='parameter']"), xpathCtx);
    if (paramObj && paramObj->nodesetval)
    {
        for (int i = 0; i < paramObj->nodesetval->nodeNr; ++i)
        {
            xmlNodePtr pNode = paramObj->nodesetval->nodeTab[i];
            std::string pname;
            xmlChar *nAttr = xmlGetProp(pNode, reinterpret_cast<const xmlChar *>("data-type:name"));
            if (!nAttr)
                nAttr = xmlGetProp(pNode, reinterpret_cast<const xmlChar *>("name"));
            if (nAttr)
            {
                pname = sanitize_workspace_string(reinterpret_cast<char *>(nAttr));
                xmlFree(nAttr);
            }

            std::string pinfix;
            xmlChar *iAttr = xmlGetProp(pNode, reinterpret_cast<const xmlChar *>("userProvidedCompInfix"));
            if (iAttr)
            {
                pinfix = sanitize_workspace_string(reinterpret_cast<char *>(iAttr));
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

    xmlXPathObjectPtr spillObj = xmlXPathNodeEval(matrixNode, reinterpret_cast<const xmlChar *>(".//*[local-name()='spillover']"), xpathCtx);
    if (spillObj && spillObj->nodesetval)
    {
        for (int i = 0; i < spillObj->nodesetval->nodeNr; ++i)
        {
            xmlNodePtr sNode = spillObj->nodesetval->nodeTab[i];
            std::string row_param;
            xmlChar *nAttr = xmlGetProp(sNode, reinterpret_cast<const xmlChar *>("data-type:parameter"));
            if (!nAttr)
                nAttr = xmlGetProp(sNode, reinterpret_cast<const xmlChar *>("parameter"));
            if (nAttr)
            {
                row_param = sanitize_workspace_string(reinterpret_cast<char *>(nAttr));
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
                xmlXPathObjectPtr coefObj = xmlXPathNodeEval(sNode, reinterpret_cast<const xmlChar *>(".//*[local-name()='coefficient']"), xpathCtx);
                if (coefObj && coefObj->nodesetval)
                {
                    for (int j = 0; j < coefObj->nodesetval->nodeNr; ++j)
                    {
                        xmlNodePtr cNode = coefObj->nodesetval->nodeTab[j];
                        std::string col_param;
                        xmlChar *cnAttr = xmlGetProp(cNode, reinterpret_cast<const xmlChar *>("data-type:parameter"));
                        if (!cnAttr)
                            cnAttr = xmlGetProp(cNode, reinterpret_cast<const xmlChar *>("parameter"));
                        if (cnAttr)
                        {
                            col_param = sanitize_workspace_string(reinterpret_cast<char *>(cnAttr));
                            xmlFree(cnAttr);
                        }

                        double val = 0.0;
                        xmlChar *vAttr = xmlGetProp(cNode, reinterpret_cast<const xmlChar *>("transforms:value"));
                        if (!vAttr)
                            vAttr = xmlGetProp(cNode, reinterpret_cast<const xmlChar *>("value"));
                        if (vAttr)
                        {
                            char *end;
                            double parsed_val = std::strtod(reinterpret_cast<char *>(vAttr), &end);
                            if (end != reinterpret_cast<char *>(vAttr)) {
                                val = parsed_val;
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

    // If still not found, prompt with native OS file dialog
    if (filename.empty() || !std::filesystem::exists(filename))
    {
#if defined(_WIN32)
        FILE *fp = _popen("powershell -NoProfile -Command \"Add-Type -AssemblyName System.Windows.Forms; $f = New-Object System.Windows.Forms.OpenFileDialog; $f.Filter = 'FlowJo Workspaces (*.wsp)|*.wsp|All Files (*.*)|*.*'; $f.Title = 'Select FlowJo Workspace'; if ($f.ShowDialog() -eq [System.Windows.Forms.DialogResult]::OK) { Write-Output $f.FileName }\" 2>NUL", "r");
        if (fp)
        {
            char path[1024] = {0};
            if (fgets(path, sizeof(path), fp))
            {
                std::string s(path);
                while (!s.empty() && (s.back() == '\n' || s.back() == '\r'))
                    s.pop_back();
                filename = s;
            }
            _pclose(fp);
        }
#elif defined(__APPLE__)
        FILE *fp = popen("osascript -e 'POSIX path of (choose file with prompt \"Select FlowJo Workspace:\" of type {\"wsp\", \"public.data\"})' 2>/dev/null", "r");
        if (fp)
        {
            char path[1024] = {0};
            if (fgets(path, sizeof(path), fp))
            {
                std::string s(path);
                while (!s.empty() && (s.back() == '\n' || s.back() == '\r'))
                    s.pop_back();
                filename = s;
            }
            pclose(fp);
        }
#endif
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
    xmlXPathObjectPtr samplesObj = xmlXPathEvalExpression(reinterpret_cast<const xmlChar *>("//SampleList/Sample"), xpathCtx);
    if (samplesObj && samplesObj->nodesetval)
    {
        for (int i = 0; i < samplesObj->nodesetval->nodeNr; ++i)
        {
            xmlNodePtr sampleNode = samplesObj->nodesetval->nodeTab[i];
            SampleData sd;

            xpathCtx->node = sampleNode;
            xmlXPathObjectPtr transObj = xmlXPathEvalExpression(reinterpret_cast<const xmlChar *>(".//*[local-name()='logicle' or local-name()='hyperlog' or local-name()='lin' or local-name()='linear' or local-name()='log' or local-name()='fasinh' or local-name()='biexp' or local-name()='biex']"), xpathCtx);
            if (transObj && transObj->nodesetval)
            {
                for (int t = 0; t < transObj->nodesetval->nodeNr; ++t)
                {
                    xmlNodePtr node = transObj->nodesetval->nodeTab[t];

                    std::string id;
                    xmlChar *idAttr = xmlGetProp(node, reinterpret_cast<const xmlChar *>("id"));
                    if (!idAttr)
                        idAttr = xmlGetProp(node, reinterpret_cast<const xmlChar *>("gating:id"));
                    if (!idAttr && node->parent)
                    {
                        idAttr = xmlGetProp(node->parent, reinterpret_cast<const xmlChar *>("id"));
                        if (!idAttr)
                            idAttr = xmlGetProp(node->parent, reinterpret_cast<const xmlChar *>("gating:id"));
                    }
                    if (idAttr)
                    {
                        id = sanitize_workspace_string(reinterpret_cast<char *>(idAttr));
                        xmlFree(idAttr);
                    }

                    if (id.empty())
                    {
                        // Fallback to checking the parameter name since they don't explicitly carry IDs
                        xmlXPathObjectPtr paramObj = xmlXPathNodeEval(node, reinterpret_cast<const xmlChar *>(".//*[local-name()='parameter']"), xpathCtx);
                        if (paramObj && paramObj->nodesetval && paramObj->nodesetval->nodeNr > 0)
                        {
                            xmlNodePtr paramNode = paramObj->nodesetval->nodeTab[0];
                            xmlChar *nameAttr = xmlGetProp(paramNode, reinterpret_cast<const xmlChar *>("data-type:name"));
                            if (!nameAttr)
                                nameAttr = xmlGetProp(paramNode, reinterpret_cast<const xmlChar *>("name"));
                            if (nameAttr)
                            {
                                id = sanitize_workspace_string(reinterpret_cast<char *>(nameAttr));
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
                        xmlChar *attr = xmlGetProp(n, reinterpret_cast<const xmlChar *>(attr_name));
                        if (attr)
                        {
                            double val = def;
                            char *end;
                            double v = std::strtod(reinterpret_cast<char *>(attr), &end);
                            if (end != reinterpret_cast<char *>(attr)) {
                                val = v;
                            }
                            xmlFree(attr);
                            return val;
                        }
                        return def;
                    };

                    std::string name = reinterpret_cast<const char *>(node->name);
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
                            double t = get_double(node, "maxRange", 262144.0);
                            double m = get_double(node, "pos", 4.418539922);
                            double a = get_double(node, "neg", 0.0);
                            double width = get_double(node, "width", -100.0);
                            double w = 0.5 * std::log10(-width);
                            transform = std::make_shared<Logicle>(t, w, m, a);
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

            xmlXPathObjectPtr sampleMatObj = xmlXPathNodeEval(sampleNode, reinterpret_cast<const xmlChar *>(".//*[local-name()='spilloverMatrix']"), xpathCtx);
            if (sampleMatObj && sampleMatObj->nodesetval && sampleMatObj->nodesetval->nodeNr > 0)
            {
                sd.spillover_matrix = parse_spillover_matrix(sampleMatObj->nodesetval->nodeTab[0], xpathCtx);
            }
            if (sampleMatObj)
                xmlXPathFreeObject(sampleMatObj);

            // Try fetching name from 'name' attribute
            xmlChar *nameAttr = xmlGetProp(sampleNode, reinterpret_cast<const xmlChar *>("name"));
            if (nameAttr)
            {
                sd.name = sanitize_workspace_string(reinterpret_cast<char *>(nameAttr));
                xmlFree(nameAttr);
            }
            else
            {
                // Fallback to finding the $FIL keyword
                xpathCtx->node = sampleNode;
                xmlXPathObjectPtr filObj = xmlXPathEvalExpression(reinterpret_cast<const xmlChar *>(".//Keyword[@name='$FIL']/@value"), xpathCtx);
                if (filObj && filObj->nodesetval && filObj->nodesetval->nodeNr > 0)
                {
                    xmlChar *content = xmlNodeGetContent(filObj->nodesetval->nodeTab[0]);
                    if (content)
                    {
                        sd.name = sanitize_workspace_string(reinterpret_cast<char *>(content));
                        xmlFree(content);
                    }
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
            xmlXPathObjectPtr varObj = xmlXPathEvalExpression(reinterpret_cast<const xmlChar *>(".//Keyword[starts-with(@name, '$P') and substring(@name, string-length(@name)) = 'N']"), xpathCtx);
            if (varObj && varObj->nodesetval)
            {
                for (int j = 0; j < varObj->nodesetval->nodeNr; ++j)
                {
                    xmlNodePtr kwNode = varObj->nodesetval->nodeTab[j];
                    xmlChar *nameAttr = xmlGetProp(kwNode, reinterpret_cast<const xmlChar *>("name"));
                    xmlChar *valAttr = xmlGetProp(kwNode, reinterpret_cast<const xmlChar *>("value"));
                    if (nameAttr && valAttr)
                    {
                        std::string varName = sanitize_workspace_string(reinterpret_cast<char *>(nameAttr));
                        sd.variables.push_back(sanitize_workspace_string(reinterpret_cast<char *>(valAttr)));

                        std::string num;
                        if (varName.size() >= 4 && varName[0] == '$' && varName[1] == 'P' && varName.back() == 'N')
                            num = varName.substr(2, varName.size() - 3);

                        const bool valid_channel = !num.empty() && std::all_of(num.begin(), num.end(), [](unsigned char c)
                                                                               { return std::isdigit(c) != 0; });
                        xmlXPathObjectPtr sObj = nullptr;
                        if (valid_channel)
                        {
                            std::string stainPath = ".//Keyword[@name='$P" + num + "S']/@value";
                            sObj = xmlXPathEvalExpression(reinterpret_cast<const xmlChar *>(stainPath.c_str()), xpathCtx);
                        }
                        if (sObj && sObj->nodesetval && sObj->nodesetval->nodeNr > 0)
                        {
                            xmlChar *sVal = xmlNodeGetContent(sObj->nodesetval->nodeTab[0]);
                            if (sVal)
                            {
                                sd.stains.push_back(sanitize_workspace_string(reinterpret_cast<char *>(sVal)));
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
            xmlXPathObjectPtr popObj = xmlXPathEvalExpression(reinterpret_cast<const xmlChar *>(".//Population"), xpathCtx);
            if (popObj && popObj->nodesetval)
            {
                for (int j = 0; j < popObj->nodesetval->nodeNr; ++j)
                {
                    xmlChar *popName = xmlGetProp(popObj->nodesetval->nodeTab[j], reinterpret_cast<const xmlChar *>("name"));
                    if (popName)
                    {
                        std::string pname = sanitize_workspace_string(reinterpret_cast<char *>(popName));
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
                                xmlXPathObjectPtr gateObj = xmlXPathEvalExpression(reinterpret_cast<const xmlChar *>(".//*[contains(local-name(), 'Gate')]"), xpathCtx);
                                if (gateObj && gateObj->nodesetval)
                                {
                                    for (int k = 0; k < gateObj->nodesetval->nodeNr; ++k)
                                    {
                                        xmlNodePtr gateNode = gateObj->nodesetval->nodeTab[k];
                                        std::string gname = reinterpret_cast<const char *>(gateNode->name);
                                        // get the graph ids
                                        if (gname == "Gate")
                                        {
                                            xmlChar *idAttr = xmlGetProp(gateNode, reinterpret_cast<const xmlChar *>("gating:id"));
                                            if (!idAttr)
                                                idAttr = xmlGetProp(gateNode, reinterpret_cast<const xmlChar *>("id"));
                                            if (idAttr)
                                            {
                                                id = reinterpret_cast<const char *>(idAttr);
                                                xmlFree(idAttr);
                                            }

                                            xmlChar *parentIdAttr = xmlGetProp(gateNode, reinterpret_cast<const xmlChar *>("gating:parent_id"));
                                            if (!parentIdAttr)
                                                parentIdAttr = xmlGetProp(gateNode, reinterpret_cast<const xmlChar *>("parent_id"));
                                                parent_id = sanitize_workspace_string(reinterpret_cast<char *>(parentIdAttr));
                                            {
                                                parent_id = reinterpret_cast<const char *>(parentIdAttr);
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
                                            xmlXPathObjectPtr vObj = xmlXPathNodeEval(gateNode, reinterpret_cast<const xmlChar *>(".//*[local-name()='vertex']"), xpathCtx);
                                            if (vObj && vObj->nodesetval)
                                            {
                                                for (int v = 0; v < vObj->nodesetval->nodeNr; ++v)
                                                {
                                                    std::vector<double> vertex;
                                                    xmlXPathObjectPtr coordObj = xmlXPathNodeEval(vObj->nodesetval->nodeTab[v], reinterpret_cast<const xmlChar *>(".//*[local-name()='coordinate']"), xpathCtx);
                                                    if (coordObj && coordObj->nodesetval)
                                                    {
                                                        for (int c = 0; c < coordObj->nodesetval->nodeNr; ++c)
                                                        {
                                                            xmlChar *valAttr = xmlGetProp(coordObj->nodesetval->nodeTab[c], reinterpret_cast<const xmlChar *>("data-type:value"));
                                                            if (!valAttr)
                                                                valAttr = xmlGetProp(coordObj->nodesetval->nodeTab[c], reinterpret_cast<const xmlChar *>("value"));
                                                            if (valAttr)
                                                            {
                                                                char *end;
                                                                double v = std::strtod(reinterpret_cast<char *>(valAttr), &end);
                                                                if (end != reinterpret_cast<char *>(valAttr)) {
                                                                    vertex.push_back(v);
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
                                            xmlXPathObjectPtr meanObj = xmlXPathNodeEval(gateNode, reinterpret_cast<const xmlChar *>(".//*[local-name()='mean']//*[local-name()='coordinate']"), xpathCtx);
                                            if (meanObj && meanObj->nodesetval)
                                            {
                                                for (int c = 0; c < meanObj->nodesetval->nodeNr; ++c)
                                                {
                                                    xmlChar *valAttr = xmlGetProp(meanObj->nodesetval->nodeTab[c], reinterpret_cast<const xmlChar *>("data-type:value"));
                                                    if (!valAttr)
                                                        valAttr = xmlGetProp(meanObj->nodesetval->nodeTab[c], reinterpret_cast<const xmlChar *>("value"));
                                                    if (valAttr)
                                                    {
                                                        char *end;
                                                        double v = std::strtod(reinterpret_cast<char *>(valAttr), &end);
                                                        if (end != reinterpret_cast<char *>(valAttr)) {
                                                            ellip->mean.push_back(v);
                                                        }
                                                        xmlFree(valAttr);
                                                    }
                                                }
                                            }
                                            if (meanObj)
                                                xmlXPathFreeObject(meanObj);
                                            xmlXPathObjectPtr rowObj = xmlXPathNodeEval(gateNode, reinterpret_cast<const xmlChar *>(".//*[local-name()='covarianceMatrix']//*[local-name()='row']"), xpathCtx);
                                            if (rowObj && rowObj->nodesetval)
                                            {
                                                for (int r = 0; r < rowObj->nodesetval->nodeNr; ++r)
                                                {
                                                    std::vector<double> row;
                                                    xmlXPathObjectPtr entryObj = xmlXPathNodeEval(rowObj->nodesetval->nodeTab[r], reinterpret_cast<const xmlChar *>(".//*[local-name()='entry']"), xpathCtx);
                                                    if (entryObj && entryObj->nodesetval)
                                                    {
                                                        for (int e = 0; e < entryObj->nodesetval->nodeNr; ++e)
                                                        {
                                                            xmlChar *valAttr = xmlGetProp(entryObj->nodesetval->nodeTab[e], reinterpret_cast<const xmlChar *>("data-type:value"));
                                                            if (!valAttr)
                                                                valAttr = xmlGetProp(entryObj->nodesetval->nodeTab[e], reinterpret_cast<const xmlChar *>("value"));
                                                            if (valAttr)
                                                            {
                                                                char *end;
                                                                double v = std::strtod(reinterpret_cast<char *>(valAttr), &end);
                                                                if (end != reinterpret_cast<char *>(valAttr)) {
                                                                    row.push_back(v);
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
                                            xmlXPathObjectPtr dimObj = xmlXPathNodeEval(gateNode, reinterpret_cast<const xmlChar *>(".//*[local-name()='dimension' or local-name()='Dimension' or local-name()='fcs-dimension']"), xpathCtx);
                                            if (dimObj && dimObj->nodesetval)
                                            {
                                                for (int d = 0; d < dimObj->nodesetval->nodeNr; ++d)
                                                {
                                                    xmlNodePtr dimNode = dimObj->nodesetval->nodeTab[d];
                                                    std::string nodeName = reinterpret_cast<const char *>(dimNode->name);
                                                    if (nodeName.find("fcs-dimension") != std::string::npos && dimNode->parent)
                                                    {
                                                        std::string parentName = reinterpret_cast<const char *>(dimNode->parent->name);
                                                        if (parentName.find("dimension") != std::string::npos || parentName.find("Dimension") != std::string::npos)
                                                        {
                                                            continue;
                                                        }
                                                    }

                                                    Gate::Dimension dim;

                                                    xmlChar *minAttr = xmlGetProp(dimNode, reinterpret_cast<const xmlChar *>("gating:min"));
                                                    if (!minAttr)
                                                        minAttr = xmlGetProp(dimNode, reinterpret_cast<const xmlChar *>("min"));
                                                    if (minAttr)
                                                    {
                                                        char *end;
                                                        double v = std::strtod(reinterpret_cast<char *>(minAttr), &end);
                                                        if (end != reinterpret_cast<char *>(minAttr)) {
                                                            dim.min_val = v;
                                                        }
                                                        xmlFree(minAttr);
                                                    }

                                                    xmlChar *maxAttr = xmlGetProp(dimNode, reinterpret_cast<const xmlChar *>("gating:max"));
                                                    if (!maxAttr)
                                                        maxAttr = xmlGetProp(dimNode, reinterpret_cast<const xmlChar *>("max"));
                                                    if (maxAttr)
                                                    {
                                                        char *end;
                                                        double v = std::strtod(reinterpret_cast<char *>(maxAttr), &end);
                                                        if (end != reinterpret_cast<char *>(maxAttr)) {
                                                            dim.max_val = v;
                                                        }
                                                        xmlFree(maxAttr);
                                                    }

                                                    xmlNodePtr paramNode = dimNode;
                                                    xmlXPathObjectPtr fcsDimObj = xmlXPathNodeEval(dimNode, reinterpret_cast<const xmlChar *>(".//*[local-name()='fcs-dimension']"), xpathCtx);
                                                    if (fcsDimObj && fcsDimObj->nodesetval && fcsDimObj->nodesetval->nodeNr > 0)
                                                    {
                                                        paramNode = fcsDimObj->nodesetval->nodeTab[0];
                                                    }
                                                    if (fcsDimObj)
                                                        xmlXPathFreeObject(fcsDimObj);

                                                    xmlChar *nAttr = xmlGetProp(paramNode, reinterpret_cast<const xmlChar *>("data-type:name"));
                                                    if (!nAttr)
                                                        nAttr = xmlGetProp(paramNode, reinterpret_cast<const xmlChar *>("name"));
                                                    if (nAttr)
                                                    {
                                                        dim.name = sanitize_workspace_string(reinterpret_cast<char *>(nAttr));
                                                        xmlFree(nAttr);
                                                    }

                                                    xmlChar *sAttr = xmlGetProp(paramNode, reinterpret_cast<const xmlChar *>("data-type:transformation-ref"));
                                                    if (!sAttr)
                                                        sAttr = xmlGetProp(paramNode, reinterpret_cast<const xmlChar *>("transformation-ref"));
                                                    if (sAttr)
                                                    {
                                                        dim.scale = sanitize_workspace_string(reinterpret_cast<char *>(sAttr));
                                                        auto it = sd.transforms.find(dim.scale);
                                                        if (it != sd.transforms.end())
                                                        {
                                                            dim.transform = it->second;
                                                        }
                                                        xmlFree(sAttr);
                                                    }

                                                    xmlChar *cAttr = xmlGetProp(paramNode, reinterpret_cast<const xmlChar *>("data-type:compensation-ref"));
                                                    if (!cAttr)
                                                        cAttr = xmlGetProp(paramNode, reinterpret_cast<const xmlChar *>("compensation-ref"));
                                                    if (cAttr)
                                                    {
                                                        dim.compensation = sanitize_workspace_string(reinterpret_cast<char *>(cAttr));
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

    const std::string safe_sample_name = sanitize_workspace_string(sample_name);
    std::string expr = "//SampleList/Sample[SampleNode/@name=" + xpath_literal(safe_sample_name + ".fcs") + " or SampleNode/@name=" + xpath_literal(safe_sample_name + ".FCS") + " or SampleNode/@name=" + xpath_literal(safe_sample_name) + " or .//Keyword[@name='$FIL' and (@value=" + xpath_literal(safe_sample_name + ".fcs") + " or @value=" + xpath_literal(safe_sample_name + ".FCS") + " or @value=" + xpath_literal(safe_sample_name) + ")]]";
    xmlXPathObjectPtr sampleObj = xmlXPathEvalExpression(reinterpret_cast<const xmlChar *>(expr.c_str()), xpathCtx);
    if (sampleObj && sampleObj->nodesetval && sampleObj->nodesetval->nodeNr > 0)
    {
        xmlNodePtr sampleNode = sampleObj->nodesetval->nodeTab[0];

        xmlXPathObjectPtr datasetObj = xmlXPathNodeEval(sampleNode, reinterpret_cast<const xmlChar *>("./DataSet"), xpathCtx);
        std::string uri;
        if (datasetObj && datasetObj->nodesetval && datasetObj->nodesetval->nodeNr > 0)
        {
            xmlChar *uriAttr = xmlGetProp(datasetObj->nodesetval->nodeTab[0], reinterpret_cast<const xmlChar *>("uri"));
            if (uriAttr)
            {
                uri = sanitize_workspace_string(reinterpret_cast<char *>(uriAttr));
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

        xmlXPathObjectPtr dpObj = xmlXPathNodeEval(sampleNode, reinterpret_cast<const xmlChar *>("./DerivedParameters"), xpathCtx);
        xmlNodePtr dpNode = nullptr;
        if (dpObj && dpObj->nodesetval && dpObj->nodesetval->nodeNr > 0)
        {
            dpNode = dpObj->nodesetval->nodeTab[0];
        }
        else
        {
            dpNode = xmlNewNode(nullptr, reinterpret_cast<const xmlChar *>("DerivedParameters"));
            xmlNodePtr keywordsNode = nullptr;
            for (xmlNodePtr child = sampleNode->children; child; child = child->next)
            {
                if (child->type == XML_ELEMENT_NODE && (xmlStrcmp(child->name, reinterpret_cast<const xmlChar *>("Keywords")) == 0 || xmlStrcmp(child->name, reinterpret_cast<const xmlChar *>("SampleNode")) == 0))
                {
                    keywordsNode = child;
                    break;
                }
            }
            if (keywordsNode)
            {
                xmlAddPrevSibling(keywordsNode, dpNode);
                xmlAddPrevSibling(keywordsNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n       ")));
            }
            else
            {
                xmlAddChild(sampleNode, dpNode);
            }
        }
        if (dpObj) xmlXPathFreeObject(dpObj);

        bool modified = false;

        xmlXPathObjectPtr pObj = xmlXPathNodeEval(dpNode, reinterpret_cast<const xmlChar *>("./DerivedParameter[@name='Laplace']"), xpathCtx);
        if (!pObj || !pObj->nodesetval || pObj->nodesetval->nodeNr == 0)
        {
            xmlNodePtr paramNode = xmlNewNode(nullptr, reinterpret_cast<const xmlChar *>("DerivedParameter"));
            xmlSetProp(paramNode, reinterpret_cast<const xmlChar *>("name"), reinterpret_cast<const xmlChar *>("Laplace"));
            xmlSetProp(paramNode, reinterpret_cast<const xmlChar *>("type"), reinterpret_cast<const xmlChar *>("importCsv"));
            xmlSetProp(paramNode, reinterpret_cast<const xmlChar *>("importFile"), reinterpret_cast<const xmlChar *>(uri.c_str()));
            xmlSetProp(paramNode, reinterpret_cast<const xmlChar *>("range"), reinterpret_cast<const xmlChar *>("1024"));
            xmlSetProp(paramNode, reinterpret_cast<const xmlChar *>("columnIndex"), reinterpret_cast<const xmlChar *>("1"));

            xmlNodePtr transNode = xmlNewNode(nullptr, reinterpret_cast<const xmlChar *>("Transform"));

            xmlNsPtr transformsNs = xmlSearchNs(doc, doc->children, reinterpret_cast<const xmlChar *>("transforms"));
            xmlNsPtr datatypeNs = xmlSearchNs(doc, doc->children, reinterpret_cast<const xmlChar *>("data-type"));

            xmlNodePtr linearNode = nullptr;
            if (transformsNs) {
                linearNode = xmlNewNode(transformsNs, reinterpret_cast<const xmlChar *>("linear"));
                xmlSetNsProp(linearNode, transformsNs, reinterpret_cast<const xmlChar *>("minRange"), reinterpret_cast<const xmlChar *>("0"));
                xmlSetNsProp(linearNode, transformsNs, reinterpret_cast<const xmlChar *>("maxRange"), reinterpret_cast<const xmlChar *>("1024"));
            } else {
                linearNode = xmlNewNode(nullptr, reinterpret_cast<const xmlChar *>("transforms:linear"));
                xmlSetProp(linearNode, reinterpret_cast<const xmlChar *>("transforms:minRange"), reinterpret_cast<const xmlChar *>("0"));
                xmlSetProp(linearNode, reinterpret_cast<const xmlChar *>("transforms:maxRange"), reinterpret_cast<const xmlChar *>("1024"));
            }
            xmlSetProp(linearNode, reinterpret_cast<const xmlChar *>("gain"), reinterpret_cast<const xmlChar *>("1"));

            xmlNodePtr pNode = nullptr;
            if (datatypeNs) {
                pNode = xmlNewNode(datatypeNs, reinterpret_cast<const xmlChar *>("parameter"));
                xmlSetNsProp(pNode, datatypeNs, reinterpret_cast<const xmlChar *>("name"), reinterpret_cast<const xmlChar *>("Laplace"));
            } else {
                pNode = xmlNewNode(nullptr, reinterpret_cast<const xmlChar *>("data-type:parameter"));
                xmlSetProp(pNode, reinterpret_cast<const xmlChar *>("data-type:name"), reinterpret_cast<const xmlChar *>("Laplace"));
            }

            xmlAddChild(linearNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n               ")));
            xmlAddChild(linearNode, pNode);
            xmlAddChild(linearNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n             ")));

            xmlAddChild(transNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n             ")));
            xmlAddChild(transNode, linearNode);
            xmlAddChild(transNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n           ")));

            xmlAddChild(paramNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n           ")));
            xmlAddChild(paramNode, transNode);
            xmlAddChild(paramNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n         ")));

            xmlAddChild(dpNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n         ")));
            xmlAddChild(dpNode, paramNode);
            xmlAddChild(dpNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n       ")));

            modified = true;
        }
        if (pObj) xmlXPathFreeObject(pObj);

        xmlXPathObjectPtr transBlockObj = xmlXPathNodeEval(sampleNode, reinterpret_cast<const xmlChar *>("./Transformations"), xpathCtx);
        xmlNodePtr transBlockNode = nullptr;
        if (transBlockObj && transBlockObj->nodesetval && transBlockObj->nodesetval->nodeNr > 0)
        {
            transBlockNode = transBlockObj->nodesetval->nodeTab[0];
        }
        else
        {
            transBlockNode = xmlNewNode(nullptr, reinterpret_cast<const xmlChar *>("Transformations"));
            xmlNodePtr keywordsNode = nullptr;
            for (xmlNodePtr child = sampleNode->children; child; child = child->next)
            {
                if (child->type == XML_ELEMENT_NODE && (xmlStrcmp(child->name, reinterpret_cast<const xmlChar *>("Keywords")) == 0 || xmlStrcmp(child->name, reinterpret_cast<const xmlChar *>("DerivedParameters")) == 0 || xmlStrcmp(child->name, reinterpret_cast<const xmlChar *>("SampleNode")) == 0))
                {
                    keywordsNode = child;
                    break;
                }
            }
            if (keywordsNode)
            {
                xmlAddPrevSibling(keywordsNode, transBlockNode);
                xmlAddPrevSibling(keywordsNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n       ")));
            }
            else
            {
                xmlAddChild(sampleNode, transBlockNode);
            }
        }
        if (transBlockObj) xmlXPathFreeObject(transBlockObj);

        std::string transExpr = "./*[local-name()='linear']/*[local-name()='parameter' and (@name='Laplace' or @*[local-name()='name']='Laplace')]";
        xmlXPathObjectPtr tCheckObj = xmlXPathNodeEval(transBlockNode, reinterpret_cast<const xmlChar *>(transExpr.c_str()), xpathCtx);
        if (!tCheckObj || !tCheckObj->nodesetval || tCheckObj->nodesetval->nodeNr == 0)
        {
            xmlNsPtr transformsNs = xmlSearchNs(doc, doc->children, reinterpret_cast<const xmlChar *>("transforms"));
            xmlNsPtr datatypeNs = xmlSearchNs(doc, doc->children, reinterpret_cast<const xmlChar *>("data-type"));

            xmlNodePtr linearNode = nullptr;
            if (transformsNs) {
                linearNode = xmlNewNode(transformsNs, reinterpret_cast<const xmlChar *>("linear"));
                xmlSetNsProp(linearNode, transformsNs, reinterpret_cast<const xmlChar *>("minRange"), reinterpret_cast<const xmlChar *>("0"));
                xmlSetNsProp(linearNode, transformsNs, reinterpret_cast<const xmlChar *>("maxRange"), reinterpret_cast<const xmlChar *>("1024"));
            } else {
                linearNode = xmlNewNode(nullptr, reinterpret_cast<const xmlChar *>("transforms:linear"));
                xmlSetProp(linearNode, reinterpret_cast<const xmlChar *>("transforms:minRange"), reinterpret_cast<const xmlChar *>("0"));
                xmlSetProp(linearNode, reinterpret_cast<const xmlChar *>("transforms:maxRange"), reinterpret_cast<const xmlChar *>("1024"));
            }
            xmlSetProp(linearNode, reinterpret_cast<const xmlChar *>("gain"), reinterpret_cast<const xmlChar *>("1"));

            xmlNodePtr pNode = nullptr;
            if (datatypeNs) {
                pNode = xmlNewNode(datatypeNs, reinterpret_cast<const xmlChar *>("parameter"));
                xmlSetNsProp(pNode, datatypeNs, reinterpret_cast<const xmlChar *>("name"), reinterpret_cast<const xmlChar *>("Laplace"));
            } else {
                pNode = xmlNewNode(nullptr, reinterpret_cast<const xmlChar *>("data-type:parameter"));
                xmlSetProp(pNode, reinterpret_cast<const xmlChar *>("data-type:name"), reinterpret_cast<const xmlChar *>("Laplace"));
            }

            xmlAddChild(linearNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n           ")));
            xmlAddChild(linearNode, pNode);
            xmlAddChild(linearNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n         ")));

            xmlAddChild(transBlockNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n         ")));
            xmlAddChild(transBlockNode, linearNode);
            xmlAddChild(transBlockNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n       ")));

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

void add_laplace_gates(const std::string &filename, const std::string &sample_name, const std::string &parent_pop_name, unsigned clusters_found, size_t laplacian_offset, const std::vector<std::vector<unsigned>>& cluster_events, const std::vector<std::string>& selected_vars)
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

    const std::string safe_sample_name = sanitize_workspace_string(sample_name);
    std::string sample_expr = "//SampleList/Sample[SampleNode/@name=" + xpath_literal(safe_sample_name + ".fcs") + " or SampleNode/@name=" + xpath_literal(safe_sample_name + ".FCS") + " or SampleNode/@name=" + xpath_literal(safe_sample_name) + " or .//Keyword[@name='$FIL' and (@value=" + xpath_literal(safe_sample_name + ".fcs") + " or @value=" + xpath_literal(safe_sample_name + ".FCS") + " or @value=" + xpath_literal(safe_sample_name) + ")]]";
    xmlXPathObjectPtr sampleObj = xmlXPathEvalExpression(reinterpret_cast<const xmlChar *>(sample_expr.c_str()), xpathCtx);
    if (!sampleObj || !sampleObj->nodesetval || sampleObj->nodesetval->nodeNr == 0) {
        if (sampleObj) xmlXPathFreeObject(sampleObj);
        xmlXPathFreeContext(xpathCtx);
        xmlFreeDoc(doc);
        return;
    }
    
    xmlNodePtr sampleNode = sampleObj->nodesetval->nodeTab[0];

    xmlNodePtr parentPopNode = nullptr;
    if (parent_pop_name == "All") {
        xmlXPathObjectPtr snObj = xmlXPathNodeEval(sampleNode, reinterpret_cast<const xmlChar *>("./SampleNode"), xpathCtx);
        if (snObj && snObj->nodesetval && snObj->nodesetval->nodeNr > 0) {
            parentPopNode = snObj->nodesetval->nodeTab[0];
        }
        if (snObj) xmlXPathFreeObject(snObj);
    } else {
        std::string pop_expr = ".//Population[@name=" + xpath_literal(sanitize_workspace_string(parent_pop_name)) + "]";
        xmlXPathObjectPtr popObj = xmlXPathNodeEval(sampleNode, reinterpret_cast<const xmlChar *>(pop_expr.c_str()), xpathCtx);
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
        if (child->type == XML_ELEMENT_NODE && xmlStrcmp(child->name, reinterpret_cast<const xmlChar *>("Subpopulations")) == 0) {
            subpopsNode = child;
        } else if (child->type == XML_ELEMENT_NODE && xmlStrcmp(child->name, reinterpret_cast<const xmlChar *>("Graph")) == 0) {
            parentGraphNode = child;
        }
    }
    if (!subpopsNode) {
        subpopsNode = xmlNewNode(nullptr, reinterpret_cast<const xmlChar *>("Subpopulations"));
        xmlAddChild(parentPopNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n           ")));
        xmlAddChild(parentPopNode, subpopsNode);
    }

    xmlNsPtr gatingNs = xmlSearchNs(doc, doc->children, reinterpret_cast<const xmlChar *>("gating"));
    xmlNsPtr datatypeNs = xmlSearchNs(doc, doc->children, reinterpret_cast<const xmlChar *>("data-type"));

    for (unsigned k = 0; k <= clusters_found + 1; ++k) {
        size_t klass = k + laplacian_offset;
        
        std::string pop_name = (k == 0) ? "Laplace.Ambiguous" : "Laplace.Cluster" + std::to_string(k);
        if (k > clusters_found)
            pop_name = "Laplace.OffScale";
        else if (cluster_events[k].size() == 0)
            continue;
        
        std::string check_expr = "./Population[@name='" + pop_name + "']";
        xmlXPathObjectPtr pCheckObj = xmlXPathNodeEval(subpopsNode, reinterpret_cast<const xmlChar *>(check_expr.c_str()), xpathCtx);
        bool exists = (pCheckObj && pCheckObj->nodesetval && pCheckObj->nodesetval->nodeNr > 0);
        if (pCheckObj) xmlXPathFreeObject(pCheckObj);
        if (exists) continue;

        xmlNodePtr popNode = xmlNewNode(nullptr, reinterpret_cast<const xmlChar *>("Population"));
        xmlSetProp(popNode, reinterpret_cast<const xmlChar *>("name"), reinterpret_cast<const xmlChar *>(pop_name.c_str()));
        xmlSetProp(popNode, reinterpret_cast<const xmlChar *>("expanded"), reinterpret_cast<const xmlChar *>("1"));
        std::string count_str = k < cluster_events.size() ? std::to_string(cluster_events[k].size()) : "0";
        xmlSetProp(popNode, reinterpret_cast<const xmlChar *>("count"), reinterpret_cast<const xmlChar *>(count_str.c_str()));

        if (parentGraphNode) {
            xmlAddChild(popNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n             ")));
            xmlNodePtr copiedGraphNode = xmlCopyNode(parentGraphNode, 1);
            xmlAddChild(popNode, copiedGraphNode);
            
            if (selected_vars.size() >= 2) {
                for (xmlNodePtr child = copiedGraphNode->children; child; child = child->next) {
                    if (child->type == XML_ELEMENT_NODE && xmlStrcmp(child->name, reinterpret_cast<const xmlChar *>("Axis")) == 0) {
                        xmlChar* dimAttr = xmlGetProp(child, reinterpret_cast<const xmlChar *>("dimension"));
                        if (dimAttr) {
                            if (xmlStrcmp(dimAttr, reinterpret_cast<const xmlChar *>("x")) == 0) {
                                const std::string safe_name = sanitize_workspace_string(selected_vars[0]);
                                xmlSetProp(child, reinterpret_cast<const xmlChar *>("name"), reinterpret_cast<const xmlChar *>(safe_name.c_str()));
                            } else if (xmlStrcmp(dimAttr, reinterpret_cast<const xmlChar *>("y")) == 0) {
                                const std::string safe_name = sanitize_workspace_string(selected_vars[1]);
                                xmlSetProp(child, reinterpret_cast<const xmlChar *>("name"), reinterpret_cast<const xmlChar *>(safe_name.c_str()));
                            }
                            xmlFree(dimAttr);
                        }
                    }
                }
            }
        }

        xmlNodePtr gateNode = xmlNewNode(nullptr, reinterpret_cast<const xmlChar *>("Gate"));
        std::string gateId = "LaplaceGate_" + sanitize_workspace_string(sample_name) + "_" + sanitize_workspace_string(parent_pop_name) + "_" + std::to_string(klass);
        xmlSetProp(gateNode, reinterpret_cast<const xmlChar *>("gating:id"), reinterpret_cast<const xmlChar *>(gateId.c_str()));

        xmlNodePtr rectGateNode = xmlNewNode(gatingNs, reinterpret_cast<const xmlChar *>("RectangleGate"));
        xmlSetProp(rectGateNode, reinterpret_cast<const xmlChar *>("eventsInside"), reinterpret_cast<const xmlChar *>("1"));
        
        xmlNodePtr dimNode = xmlNewNode(gatingNs, reinterpret_cast<const xmlChar *>("dimension"));
        std::string min_val = std::to_string(1 + 4 * klass);
        std::string max_val = std::to_string(3 + 4 * klass);
        xmlSetProp(dimNode, reinterpret_cast<const xmlChar *>("gating:min"), reinterpret_cast<const xmlChar *>(min_val.c_str()));
        xmlSetProp(dimNode, reinterpret_cast<const xmlChar *>("gating:max"), reinterpret_cast<const xmlChar *>(max_val.c_str()));

        xmlNodePtr fcsDimNode = xmlNewNode(datatypeNs, reinterpret_cast<const xmlChar *>("fcs-dimension"));
        xmlSetProp(fcsDimNode, reinterpret_cast<const xmlChar *>("data-type:name"), reinterpret_cast<const xmlChar *>("Laplace"));

        xmlAddChild(dimNode, fcsDimNode);
        xmlAddChild(rectGateNode, dimNode);
        xmlAddChild(gateNode, rectGateNode);
        xmlAddChild(popNode, gateNode);
        
        xmlAddChild(subpopsNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n             ")));
        xmlAddChild(subpopsNode, popNode);
    }
    xmlAddChild(subpopsNode, xmlNewText(reinterpret_cast<const xmlChar *>("\n           ")));

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

    state.smoothing_str = "0.01";
    state.threshold_str = "0.001";
    state.max_clusters_str = "12";
    state.min_events_str = "0";
    state.kld_norm_str = "0.04";
    state.kld_exp_str = "0.2";
    state.min_cluster_rel_str = "0.0";
    state.tolerance_str = "0.01";

    int analysis_choice = 1;
    std::string settings_path = get_settings_path();
    if (!settings_path.empty() && std::filesystem::exists(settings_path)) {
        try {
            std::ifstream in(settings_path);
            json j;
            in >> j;
            state.smoothing_str = j.value("smoothing", "0.01");
            state.threshold_str = j.value("threshold", "0.001");
            state.max_clusters_str = j.value("max_clusters", "12");
            state.min_events_str = j.value("min_events", "0");
            state.kld_norm_str = j.value("kld_norm", "0.04");
            state.kld_exp_str = j.value("kld_exp", "0.2");
            state.min_cluster_rel_str = j.value("min_cluster_rel", "0.0");
            state.tolerance_str = j.value("tolerance", "0.01");
            analysis_choice = j.value("analysis_choice", 1);
        } catch (...) {}
    }

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
    auto input_kld_norm = Input(&state.kld_norm_str, "0.04");
    auto input_kld_exp = Input(&state.kld_exp_str, "0.2");

    auto input_max_clusters = Input(&state.max_clusters_str, "12");
    auto input_min_events = Input(&state.min_events_str, "0");
    auto input_min_cluster_rel = Input(&state.min_cluster_rel_str, "0.0");
    auto input_tolerance = Input(&state.tolerance_str, "0.01");

    auto parameters_container = Container::Vertical({
        input_smoothing,
        input_min_events,
        input_min_cluster_rel
    });

    auto laplace_container = Container::Vertical({
        input_threshold,
        input_max_clusters
    });

    auto klds_container = Container::Vertical({
        input_tolerance,
        input_kld_norm,
        input_kld_exp
    });

    auto main_layout = Container::Horizontal({sample_container,
                                              Maybe(var_container, [&]
                                                    { return num_samples_selected > 0; }),
                                              Maybe(pop_container, [&]
                                                    { return num_vars_selected >= 2; })});

    auto bottom_container = Container::Horizontal({choice_handled,
                                                   parameters_container,
                                                   laplace_container,
                                                   klds_container});

    auto top_level = Container::Vertical({main_layout,
                                          bottom_container});

    auto screen = ScreenInteractive::TerminalOutput();
    auto top_level_handled = CatchEvent(top_level, [&, main_layout, bottom_container](ftxui::Event e)
                                        {
        if (e == ftxui::Event::Return) {
            screen.ExitLoopClosure()();
            return true;
        }
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

        auto parameters_win = window(text(" Parameters "), vbox({
            hbox({text("Smoothing: "), input_smoothing->Render() | size(WIDTH, EQUAL, 6)}),
            hbox({text("Min Abs:   "), input_min_events->Render() | size(WIDTH, EQUAL, 6)}),
            hbox({text("Min Rel:   "), input_min_cluster_rel->Render() | size(WIDTH, EQUAL, 6)})
        }));

        auto laplace_win = window(text(" Laplace "), vbox({
            hbox({text("Threshold:    "), input_threshold->Render() | size(WIDTH, EQUAL, 6)}),
            hbox({text("Max Clusters: "), input_max_clusters->Render() | size(WIDTH, EQUAL, 6)})
        }));

        auto epp_win = window(text(" EPP "), vbox({
            hbox({text("Tolerance: "), input_tolerance->Render() | size(WIDTH, EQUAL, 6)}),
            hbox({text("KLD Norm:  "), input_kld_norm->Render() | size(WIDTH, EQUAL, 6)}),
            hbox({text("KLD Exp:   "), input_kld_exp->Render() | size(WIDTH, EQUAL, 6)})
        }));
            
        auto pop_win = num_vars_selected >= 2 
            ? window(text(" Populations "), pop_container->Render() | vscroll_indicator | frame) | size(WIDTH, EQUAL, 32)
            : emptyElement();
            
        return vbox({
            text(" Leonard - " + sanitize_workspace_string(std::filesystem::path(ws.filename).stem().string()) + " ") | bold | hcenter,
            separator(),
            hbox({ sample_win, var_win, pop_win }) | size(HEIGHT, LESS_THAN, 20),
            separator(),
            hbox({
                window(text(" Analysis Method "), choice_handled->Render()) | flex,
                parameters_win,
                laplace_win,
                epp_win
            }),
            separator(),
            text(" Space: Select | Arrows: Navigate | Tab: Switch Section | Enter: Confirm | Esc: Cancel ") | hcenter
        }) | border; });

    {
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

        switch (state.variables.size())
        {
        case 3:
            state.grid_size = 128;
            break;
        case 4:
            state.grid_size = 48;
            break;
        default:
            state.grid_size = 256;
            break;
        }

        char* end;
        float sf = std::strtof(state.smoothing_str.c_str(), &end);
        if (end != state.smoothing_str.c_str()) state.smoothing = sf;
        else state.smoothing = 0.01f;

        float tf = std::strtof(state.threshold_str.c_str(), &end);
        if (end != state.threshold_str.c_str()) state.threshold = tf;
        else state.threshold = 0.001f;

        unsigned long mc = std::strtoul(state.max_clusters_str.c_str(), &end, 10);
        if (end != state.max_clusters_str.c_str()) state.max_clusters = mc;
        else state.max_clusters = 12;

        unsigned long long me = std::strtoull(state.min_events_str.c_str(), &end, 10);
        if (end != state.min_events_str.c_str()) state.min_events = me;
        else state.min_events = 0;

        float kn = std::strtof(state.kld_norm_str.c_str(), &end);
        if (end != state.kld_norm_str.c_str()) state.kld_norm = kn;
        else state.kld_norm = 0.04f;

        float ke = std::strtof(state.kld_exp_str.c_str(), &end);
        if (end != state.kld_exp_str.c_str()) state.kld_exp = ke;
        else state.kld_exp = 0.2f;

        float mcr = std::strtof(state.min_cluster_rel_str.c_str(), &end);
        if (end != state.min_cluster_rel_str.c_str()) state.min_cluster_rel = mcr;
        else state.min_cluster_rel = 0.0f;

        float tol = std::strtof(state.tolerance_str.c_str(), &end);
        if (end != state.tolerance_str.c_str()) state.tolerance = tol;
        else state.tolerance = 0.01f;

        if (!settings_path.empty()) {
            json j;
            j["smoothing"] = state.smoothing_str;
            j["threshold"] = state.threshold_str;
            j["max_clusters"] = state.max_clusters_str;
            j["min_events"] = state.min_events_str;
            j["kld_norm"] = state.kld_norm_str;
            j["kld_exp"] = state.kld_exp_str;
            j["min_cluster_rel"] = state.min_cluster_rel_str;
            j["tolerance"] = state.tolerance_str;
            j["analysis_choice"] = state.analysis_choice;
            std::ofstream out(settings_path);
            out << j.dump(4);
        }
    }

    return state;
}