#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <fstream>
#include <sstream>
#include <map>
#include <cmath>
#include <algorithm>
#include <libxml/parser.h>
#include <libxslt/xslt.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>
#include <filesystem>

#include "Samples.hpp"

Variable::Variable() : scale(Scale::unknown) {}

Variable::~Variable() = default;

void Variable::evaluate() {}

double Variable::transform_value(double value) const
{
    if (transform)
        return transform->scale(value);
    return value;
}

double Variable::inverse_transform(double value) const
{
    if (transform)
        return transform->inverse(value);
    return value;
}

void DataSet::for_each_class(const Classification &cls, std::function<void(unsigned short)> func)
{
    if (cls.classifications)
    {
        for (auto val : *cls.classifications)
        {
            func(val);
        }
    }
}

void DataSet::for_each_member(const Subset &sub, std::function<void(bool)> func)
{
    if (sub.membership)
    {
        for (auto val : *sub.membership)
        {
            func(val);
        }
    }
}

void DataSet::setup_float_variables(const std::vector<std::string> &headers)
{
    variables.clear();
    variable.clear();
    subset.clear();
    classification.clear();
    for (const auto &h : headers)
    {
        auto &var = variable[h];
        var.data = std::make_shared<std::vector<float>>();
        variables.push_back(&var);
    }
}

bool DataSet::read_csv(const std::string &filename, char delimiter)
{
    std::ifstream in(filename);
    if (!in.is_open())
        return false;

    std::string line;
    if (!std::getline(in, line))
        return false;

    // Parse the header for variable names
    std::vector<std::string> headers;
    std::stringstream ss(line);
    std::string token;
    while (std::getline(ss, token, delimiter))
    {
        token.erase(0, token.find_first_not_of(" \r\n\t"));
        token.erase(token.find_last_not_of(" \r\n\t") + 1);
        headers.push_back(token);
    }

    setup_float_variables(headers);
    size_t num_cols = headers.size();

    // Parse data rows
    while (std::getline(in, line))
    {
        if (line.empty() || line.find_first_not_of(" \r\n\t") == std::string::npos)
            continue;

        std::stringstream line_ss(line);
        for (size_t i = 0; i < num_cols; ++i)
        {
            if (!std::getline(line_ss, token, delimiter))
                break;
            try
            {
                float val = std::stof(token);
                variables[i]->data->push_back(val);
            }
            catch (...)
            {
                variables[i]->data->push_back(std::nanf(""));
            }
        }
        ++_size;
    }
    return true;
}

bool DataSet::read_xml_xslt(const std::string &xml_file, const std::string &xsl_file, char delimiter)
{
    xmlDocPtr xml_doc = xmlParseFile(xml_file.c_str());
    if (!xml_doc)
        return false;

    xsltStylesheetPtr xsl_doc = xsltParseStylesheetFile((const xmlChar *)xsl_file.c_str());
    if (!xsl_doc)
    {
        xmlFreeDoc(xml_doc);
        return false;
    }

    xmlDocPtr res_doc = xsltApplyStylesheet(xsl_doc, xml_doc, NULL);
    std::string tmp_csv = xml_file + ".tmp.csv";
    bool success = false;

    if (res_doc)
    {
        if (xsltSaveResultToFilename(tmp_csv.c_str(), res_doc, xsl_doc, 0) != -1)
        {
            success = true;
        }
        xmlFreeDoc(res_doc);
    }
    xsltFreeStylesheet(xsl_doc);
    xmlFreeDoc(xml_doc);

    // Parse the CSV output and cleanup
    if (success)
    {
        success = read_csv(tmp_csv, delimiter);
        std::remove(tmp_csv.c_str());
    }
    return success;
}

bool DataSet::read_binary(const std::string &filename, const std::vector<std::string> &headers)
{
    std::ifstream in(filename, std::ios::binary | std::ios::ate);
    if (!in.is_open())
        return false;

    setup_float_variables(headers);
    size_t num_cols = headers.size();
    if (num_cols == 0)
        return false;

    size_t file_size = in.tellg();
    in.seekg(0, std::ios::beg);

    size_t num_floats = file_size / sizeof(float);
    size_t num_rows = num_floats / num_cols;

    std::vector<float> row(num_cols);
    for (size_t r = 0; r < num_rows; ++r)
    {
        in.read(reinterpret_cast<char *>(row.data()), num_cols * sizeof(float));
        for (size_t c = 0; c < num_cols; ++c)
        {
            variables[c]->data->push_back(row[c]);
        }
    }
    _size = num_rows;
    return in.good();
}

bool DataSet::read_fcs(const std::string &filename)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in.is_open())
        return false;

    char version[11] = {0};
    in.read(version, 10); // Check FCS2.0, FCS3.0, etc.

    auto read_offset = [&in]() -> size_t
    {
        char buf[9] = {0};
        in.read(buf, 8);
        std::string s(buf);
        s.erase(0, s.find_first_not_of(' '));
        return s.empty() ? 0 : std::stoull(s);
    };

    size_t text_start = read_offset();
    size_t text_end = read_offset();
    size_t data_start = read_offset();
    size_t data_end = read_offset();

    // Parse the TEXT segment for basic structure dimensions and labels
    in.seekg(text_start, std::ios::beg);
    size_t text_size = text_end - text_start + 1;
    std::vector<char> text_seg(text_size);
    in.read(text_seg.data(), text_size);

    char delimiter = text_seg[0];
    std::map<std::string, std::string> metadata;
    std::string key, value;
    bool is_key = true;

    for (size_t i = 1; i < text_size; ++i)
    {
        if (text_seg[i] == delimiter)
        {
            if (i + 1 < text_size && text_seg[i + 1] == delimiter)
            { // Escaped delimiter
                (is_key ? key : value) += delimiter;
                i++;
            }
            else
            {
                if (is_key)
                    is_key = false;
                else
                {
                    std::transform(key.begin(), key.end(), key.begin(), ::toupper);

                    size_t start = value.find_first_not_of(" \r\n\t");
                    if (start == std::string::npos)
                        value.clear();
                    else
                    {
                        value.erase(0, start);
                        value.erase(value.find_last_not_of(" \r\n\t") + 1);
                    }

                    metadata[key] = value;
                    key.clear();
                    value.clear();
                    is_key = true;
                }
            }
        }
        else
        {
            (is_key ? key : value) += text_seg[i];
        }
    }
    if (!key.empty())
    {
        std::transform(key.begin(), key.end(), key.begin(), ::toupper);

        size_t start = value.find_first_not_of(" \r\n\t");
        if (start == std::string::npos)
            value.clear();
        else
        {
            value.erase(0, start);
            value.erase(value.find_last_not_of(" \r\n\t") + 1);
        }

        metadata[key] = value;
    }

    // Fallback checks for > 99,999,999 byte offsets in FCS 3.0+
    if (data_start == 0 && metadata.count("$BEGINDATA"))
        data_start = std::stoull(metadata["$BEGINDATA"]);
    if (data_end == 0 && metadata.count("$ENDDATA"))
        data_end = std::stoull(metadata["$ENDDATA"]);

    int num_params = metadata.count("$PAR") ? std::stoi(metadata["$PAR"]) : 0;
    size_t tot_events = metadata.count("$TOT") ? std::stoull(metadata["$TOT"]) : 0;

    std::vector<std::string> headers;
    for (int i = 1; i <= num_params; ++i)
    {
        std::string p_name = "$P" + std::to_string(i) + "N";
        std::string p_stain = "$P" + std::to_string(i) + "S";
        // needed because spillover matrix does it
        std::string mangle = metadata[p_name];
        std::replace(mangle.begin(), mangle.end(), '/', '_');
        headers.push_back(mangle);
    }

    setup_float_variables(headers);

    // Check if little or big endian. Assuming host CPU is little-endian x86_64
    bool swap_bytes = (metadata.count("$BYTEORD") && metadata["$BYTEORD"] == "4,3,2,1");

    // Extract raw floats from DATA segment
    in.seekg(data_start, std::ios::beg);
    for (size_t e = 0; e < tot_events; ++e)
    {
        for (int p = 0; p < num_params; ++p)
        {
            float val = 0.0f;
            in.read(reinterpret_cast<char *>(&val), sizeof(float));
            if (swap_bytes)
            {
                char *v_bytes = reinterpret_cast<char *>(&val);
                std::swap(v_bytes[0], v_bytes[3]);
                std::swap(v_bytes[1], v_bytes[2]);
            }
            variables[p]->data->push_back(val);
        }
    }
    _size = tot_events;

    // Check for auxiliary .clr, .csv, and .truth files associated with this .fcs file
    auto read_aux = [this, tot_events](const std::filesystem::path &path, bool is_bool)
    {
        std::ifstream in(path);
        if (!in.is_open())
            return;

        std::string line;
        if (!std::getline(in, line))
            return;

        std::vector<std::string> aux_headers;
        std::stringstream ss(line);
        std::string token;
        while (std::getline(ss, token, ','))
        {
            token.erase(0, token.find_first_not_of(" \r\n\t"));
            token.erase(token.find_last_not_of(" \r\n\t") + 1);
            aux_headers.push_back(token);
        }

        size_t num_cols = aux_headers.size();
        std::vector<std::shared_ptr<std::vector<bool>>> bool_cols;
        std::vector<std::shared_ptr<std::vector<unsigned short>>> class_cols;

        for (const auto &h : aux_headers)
        {
            if (is_bool)
            {
                auto &sub = this->subset[h];
                sub.membership = std::make_shared<std::vector<bool>>();
                sub.membership->reserve(tot_events);
                bool_cols.push_back(sub.membership);
            }
            else
            {
                auto &cls = this->classification[h];
                cls.classifications = std::make_shared<std::vector<unsigned short>>();
                cls.classifications->reserve(tot_events);
                class_cols.push_back(cls.classifications);
            }
        }

        while (std::getline(in, line))
        {
            if (line.empty() || line.find_first_not_of(" \r\n\t") == std::string::npos)
                continue;
            std::stringstream line_ss(line);
            for (size_t i = 0; i < num_cols; ++i)
            {
                std::string t;
                if (!std::getline(line_ss, t, ','))
                    break;
                t.erase(0, t.find_first_not_of(" \r\n\t"));
                t.erase(t.find_last_not_of(" \r\n\t") + 1);

                if (is_bool)
                {
                    bool val = false;
                    if (!t.empty())
                    {
                        try
                        {
                            val = (std::stoi(t) != 0);
                        }
                        catch (...)
                        {
                        }
                    }
                    bool_cols[i]->push_back(val);
                }
                else
                {
                    unsigned short val = 0;
                    if (!t.empty())
                    {
                        try
                        {
                            val = static_cast<unsigned short>(std::stoul(t));
                        }
                        catch (...)
                        {
                        }
                    }
                    class_cols[i]->push_back(val);
                }
            }
        }
    };

    std::filesystem::path fcs_path(filename);

    std::filesystem::path clr_path = fcs_path;
    clr_path.replace_extension(".clr");
    if (std::filesystem::exists(clr_path))
    {
        read_aux(clr_path, true);
    }

    std::filesystem::path csv_path = fcs_path;
    csv_path.replace_extension(".csv");
    if (std::filesystem::exists(csv_path))
    {
        read_aux(csv_path, false);
    }

    std::filesystem::path truth_path = fcs_path;
    truth_path.replace_extension(".truth");
    if (std::filesystem::exists(truth_path))
    {
        read_aux(truth_path, false);
    }

    std::filesystem::path len_path = fcs_path;
    len_path.replace_extension(".len");
    if (std::filesystem::exists(len_path))
    {
        read_aux(len_path, false);
    }

    return true;
}

bool DataSet::read(const std::string &filename, const std::vector<std::string> &headers)
{
    auto do_read = [this, &headers](const std::string &fname, const std::string &ext) -> bool
    {
        std::string e = ext;
        std::transform(e.begin(), e.end(), e.begin(), ::tolower);
        if (e == ".fcs")
            return read_fcs(fname);
        if (e == ".csv")
            return read_csv(fname, ',');
        if (e == ".tsv" || e == ".txt")
            return read_csv(fname, '\t');
        if (e == ".bin" || e == ".dat")
            return read_binary(fname, headers);
        return false;
    };

    std::filesystem::path p(filename);
    if (p.has_extension())
    {
        std::string ext = p.extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext == ".fcs" || ext == ".csv" || ext == ".tsv" || ext == ".txt" || ext == ".bin" || ext == ".dat")
        {
            return do_read(filename, ext);
        }
    }

    const char *exts[] = {".fcs", ".csv", ".tsv", ".txt", ".bin", ".dat"};
    for (const char *ext : exts)
    {
        if (std::filesystem::exists(filename + ext))
        {
            if (do_read(filename + ext, ext))
                return true;
        }
    }

    return false;
}

size_t DataSet::size() const
{
    return _size;
}
