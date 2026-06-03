#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cstdlib>
#include <filesystem>
#include <algorithm>
#include <cctype>

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

#include "Samples.hpp"

int main(int argc, char* argv[]) {
    std::string filename;
    if (argc >= 2) {
        filename = argv[1];
        if (!std::filesystem::exists(filename)) {
            if (std::filesystem::exists(filename + ".wsp")) filename += ".wsp";
        }
    }
    
    // If not provided or missing, locate the most recent .wsp file
    if (filename.empty() || !std::filesystem::exists(filename)) {
        std::filesystem::file_time_type latest_time = std::filesystem::file_time_type::min();
        for (const auto& entry : std::filesystem::directory_iterator(".")) {
            if (entry.is_regular_file() && entry.path().extension() == ".wsp") {
                auto mtime = entry.last_write_time();
                if (mtime > latest_time || filename.empty()) {
                    latest_time = mtime;
                    filename = entry.path().string();
                }
            }
        }
    }

    if (filename.empty() || !std::filesystem::exists(filename)) {
        std::cerr << "Error: No workspace (.wsp) file found or specified.\n";
        return 1;
    }

    // --- XML Parsing ---
    xmlDocPtr doc = xmlParseFile(filename.c_str());
    if (!doc) {
        std::cerr << "Error: Could not parse workspace file: " << filename << "\n";
        return 1;
    }
    xmlXPathContextPtr xpathCtx = xmlXPathNewContext(doc);
    if (!xpathCtx) {
        std::cerr << "Error: Could not create XPath context.\n";
        xmlFreeDoc(doc);
        return 1;
    }

    struct SampleData {
        std::string name;
        std::vector<std::string> variables;
        std::vector<std::string> stains;
        std::vector<std::string> populations;
        bool selected = false;
        bool enabled = true;
    };
    
    std::vector<SampleData> samples;
    std::vector<std::string> all_populations;

    // Fetch Samples or SampleNodes
    xmlXPathObjectPtr samplesObj = xmlXPathEvalExpression((const xmlChar*)"//SampleList/Sample", xpathCtx);
    if (samplesObj && samplesObj->nodesetval) {
        for (int i = 0; i < samplesObj->nodesetval->nodeNr; ++i) {
            xmlNodePtr sampleNode = samplesObj->nodesetval->nodeTab[i];
            SampleData sd;
            
            // Try fetching name from 'name' attribute
            xmlChar* nameAttr = xmlGetProp(sampleNode, (const xmlChar*)"name");
            if (nameAttr) {
                sd.name = (char*)nameAttr;
                xmlFree(nameAttr);
            } else {
                // Fallback to finding the $FIL keyword
                xpathCtx->node = sampleNode;
                xmlXPathObjectPtr filObj = xmlXPathEvalExpression((const xmlChar*)".//Keyword[@name='$FIL']/@value", xpathCtx);
                if (filObj && filObj->nodesetval && filObj->nodesetval->nodeNr > 0) {
                    sd.name = (char*)xmlNodeGetContent(filObj->nodesetval->nodeTab[0]);
                } else {
                    sd.name = "Sample_" + std::to_string(i);
                }
                if (filObj) xmlXPathFreeObject(filObj);
            }
            
            size_t fcs_pos = sd.name.find(".fcs");
            if (fcs_pos != std::string::npos) sd.name.erase(fcs_pos, 4);
            size_t FCS_pos = sd.name.find(".FCS");
            if (FCS_pos != std::string::npos) sd.name.erase(FCS_pos, 4);
            
            // Collect $PnN keywords (variables) and $PnS keywords (stains)
            xpathCtx->node = sampleNode;
            xmlXPathObjectPtr varObj = xmlXPathEvalExpression((const xmlChar*)".//Keyword[starts-with(@name, '$P') and substring(@name, string-length(@name)) = 'N']", xpathCtx);
            if (varObj && varObj->nodesetval) {
                for (int j = 0; j < varObj->nodesetval->nodeNr; ++j) {
                    xmlNodePtr kwNode = varObj->nodesetval->nodeTab[j];
                    xmlChar* nameAttr = xmlGetProp(kwNode, (const xmlChar*)"name");
                    xmlChar* valAttr = xmlGetProp(kwNode, (const xmlChar*)"value");
                    if (nameAttr && valAttr) {
                        std::string varName = (char*)nameAttr;
                        sd.variables.push_back((char*)valAttr);
                        
                        std::string num = varName.substr(2, varName.length() - 3);
                        std::string stainPath = ".//Keyword[@name='$P" + num + "S']/@value";
                        xmlXPathObjectPtr sObj = xmlXPathEvalExpression((const xmlChar*)stainPath.c_str(), xpathCtx);
                        if (sObj && sObj->nodesetval && sObj->nodesetval->nodeNr > 0) {
                            xmlChar* sVal = xmlNodeGetContent(sObj->nodesetval->nodeTab[0]);
                            if (sVal) {
                                sd.stains.push_back((char*)sVal);
                                xmlFree(sVal);
                            } else {
                                sd.stains.push_back("");
                            }
                        } else {
                            sd.stains.push_back("");
                        }
                        if (sObj) xmlXPathFreeObject(sObj);
                    }
                    if (nameAttr) xmlFree(nameAttr);
                    if (valAttr) xmlFree(valAttr);
                }
            }
            if (varObj) xmlXPathFreeObject(varObj);
            
            // Collect Populations
            xmlXPathObjectPtr popObj = xmlXPathEvalExpression((const xmlChar*)".//Population", xpathCtx);
            if (popObj && popObj->nodesetval) {
                for (int j = 0; j < popObj->nodesetval->nodeNr; ++j) {
                    xmlChar* popName = xmlGetProp(popObj->nodesetval->nodeTab[j], (const xmlChar*)"name");
                    if (popName) {
                        std::string pname = (char*)popName;
                        pname.erase(std::remove_if(pname.begin(), pname.end(), [](unsigned char c) {
                            return c < 32 || c == 127;
                        }), pname.end());
                        auto start = pname.find_first_not_of(" ");
                        if (start != std::string::npos) {
                            pname = pname.substr(start);
                        } else {
                            pname = "";
                        }
                        auto end = pname.find_last_not_of(" ");
                        if (end != std::string::npos) {
                            pname = pname.substr(0, end + 1);
                        }

                        if (!pname.empty()) {
                            if (std::find(sd.populations.begin(), sd.populations.end(), pname) == sd.populations.end()) {
                                sd.populations.push_back(pname);
                            }
                            if (std::find(all_populations.begin(), all_populations.end(), pname) == all_populations.end()) {
                                all_populations.push_back(pname);
                            }
                        }
                        xmlFree(popName);
                    }
                }
            }
            if (popObj) xmlXPathFreeObject(popObj);
            
            samples.push_back(sd);
        }
    }
    if (samplesObj) xmlXPathFreeObject(samplesObj);
    xmlXPathFreeContext(xpathCtx);
    xmlFreeDoc(doc);
    xmlCleanupParser();
    
    std::sort(all_populations.begin(), all_populations.end());

    // Consolidate unique variables across all samples
    std::vector<std::string> all_variables;
    for (const auto& s : samples) {
        for (const auto& v : s.variables) {
            if (std::find(all_variables.begin(), all_variables.end(), v) == all_variables.end()) {
                all_variables.push_back(v);
            }
        }
    }
    std::sort(all_variables.begin(), all_variables.end());

    // --- FTXUI Interface ---
    using namespace ftxui;

    auto var_states = std::make_unique<bool[]>(all_variables.size());
    for (size_t i = 0; i < all_variables.size(); ++i) var_states[i] = false;
    
    auto pop_states = std::make_unique<bool[]>(all_populations.size());
    for (size_t i = 0; i < all_populations.size(); ++i) pop_states[i] = false;
    
    int num_samples_selected = 0;
    SampleData* selected_sample = nullptr;
    int num_vars_selected = 0;
    int num_pops_selected = 0;
    std::vector<std::string> selected_pops;

    auto sample_container = Container::Vertical({});
    for (size_t i = 0; i < samples.size(); ++i) {
        auto cb = Checkbox(&samples[i].name, &samples[i].selected);
        auto cb_styled = Renderer(cb, [cb, &samples, i] {
            auto el = cb->Render();
            if (!samples[i].enabled) el = el | color(Color::GrayDark);
            return el;
        });
        auto cb_handled = CatchEvent(cb_styled, [&samples, i](Event e) {
            if (!samples[i].enabled) {
                if (e == Event::Character(' ') || e == Event::Return) return true; // block interaction
            }
            return false;
        });
        sample_container->Add(cb_handled);
    }

    auto var_container = Container::Vertical({});
    for (size_t i = 0; i < all_variables.size(); ++i) {
        auto cb = Checkbox(&all_variables[i], &var_states[i]);
        var_container->Add(Maybe(cb, [&, i] {
            if (num_samples_selected == 0 || !selected_sample) return false;
            return std::find(selected_sample->variables.begin(), selected_sample->variables.end(), all_variables[i]) != selected_sample->variables.end();
        }));
    }

    auto pop_container = Container::Vertical({});
    for (size_t i = 0; i < all_populations.size(); ++i) {
        auto cb = Checkbox(&all_populations[i], &pop_states[i]);
        pop_container->Add(Maybe(cb, [&, i] {
            if (num_samples_selected == 0 || !selected_sample) return false;
            return std::find(selected_sample->populations.begin(), selected_sample->populations.end(), all_populations[i]) != selected_sample->populations.end();
        }));
    }

    int analysis_choice = 0;
    std::vector<std::string> choices = {
        "Complex sample with many dimensions, analyze with two dimensional methods that can be sorted and visualized.",
        "Sample that cannot be resolved in two dimensions, try higher dimensional methods that cannot be sorted or visualized."
    };
    auto choice_container = Radiobox(&choices, &analysis_choice);
    auto choice_handled = CatchEvent(choice_container, [&](Event e) {
        if (num_vars_selected > 4 && analysis_choice == 0) {
            if (e == Event::ArrowDown || e == Event::Character('j')) return true; // block switching to 2nd option
        }
        return false;
    });

    auto main_layout = Container::Horizontal({
        sample_container,
        Maybe(var_container, [&] { return num_samples_selected > 0; }),
        Maybe(pop_container, [&] { return num_vars_selected >= 2; })
    });

    auto top_level = Container::Vertical({
        main_layout,
        choice_handled
    });

    auto top_level_handled = CatchEvent(top_level, [&, main_layout, choice_handled](Event e) {
        if (e == Event::Tab) {
            if (main_layout->Focused()) {
                choice_handled->TakeFocus();
            } else {
                main_layout->TakeFocus();
            }
            return true;
        }
        return false;
    });

    auto renderer = Renderer(top_level_handled, [&] {
        // Pre-render state resolution 
        num_samples_selected = 0;
        selected_sample = nullptr;
        for (auto& s : samples) {
            if (s.selected) {
                num_samples_selected++;
                selected_sample = &s;
            }
        }
        
        num_vars_selected = 0;
        for (size_t i = 0; i < all_variables.size(); ++i) {
            if (var_states[i] && selected_sample && 
                std::find(selected_sample->variables.begin(), selected_sample->variables.end(), all_variables[i]) != selected_sample->variables.end()) {
                num_vars_selected++;
            } else {
                var_states[i] = false; // ensure variables hidden by a new sample are deselected
            }
        }
        
        num_pops_selected = 0;
        selected_pops.clear();
        for (size_t i = 0; i < all_populations.size(); ++i) {
            bool in_sample = selected_sample && std::find(selected_sample->populations.begin(), selected_sample->populations.end(), all_populations[i]) != selected_sample->populations.end();
            if (pop_states[i] && in_sample) {
                num_pops_selected++;
                selected_pops.push_back(all_populations[i]);
            } else {
                pop_states[i] = false;
            }
        }
        
        if (num_vars_selected > 4) analysis_choice = 0;
        
        for (auto& s : samples) {
            if (num_pops_selected > 0) {
                bool has_all = true;
                for (const auto& p : selected_pops) {
                    if (std::find(s.populations.begin(), s.populations.end(), p) == s.populations.end()) {
                        has_all = false; break;
                    }
                }
                s.enabled = has_all;
            } else if (num_samples_selected > 0) {
                s.enabled = s.selected;
            } else {
                s.enabled = true;
            }
        }

        Elements stain_elements;
        if (selected_sample) {
            for (size_t i = 0; i < all_variables.size(); ++i) {
                auto it = std::find(selected_sample->variables.begin(), selected_sample->variables.end(), all_variables[i]);
                if (it != selected_sample->variables.end()) {
                    int idx = std::distance(selected_sample->variables.begin(), it);
                    stain_elements.push_back(text(selected_sample->stains[idx]) | size(HEIGHT, EQUAL, 1));
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
            
        auto pop_win = num_vars_selected >= 2 
            ? window(text(" Populations "), pop_container->Render() | vscroll_indicator | frame) | size(WIDTH, EQUAL, 30)
            : emptyElement();
            
        return vbox({
            text(" Leonard - " + std::filesystem::path(filename).stem().string() + " ") | bold | hcenter,
            separator(),
            hbox({ sample_win, var_win, pop_win }) | size(HEIGHT, LESS_THAN, 20),
            separator(),
            window(text(" Analysis Method "), choice_handled->Render()),
            separator(),
            text(" Space: Select | Arrows: Navigate | Tab: Switch Section | Enter: Confirm | Esc: Cancel ") | hcenter
        }) | border;
    });

    bool is_cancelled = false;
    {
        auto screen = ScreenInteractive::TerminalOutput();
        auto main_container = CatchEvent(renderer, [&](ftxui::Event event) {
            if (event == ftxui::Event::Return) {
                screen.ExitLoopClosure()();
                return true;
            }
            if (event == ftxui::Event::Escape) {
                is_cancelled = true;
                screen.ExitLoopClosure()();
                return true;
            }
            return false;
        });

        screen.Loop(main_container);
    }

    std::cout << "\x1B[2J\x1B[H";
    std::cout.flush();

#ifndef _WIN32
    // Flush any stray terminal responses (e.g., cursor position) received during the Escape timeout
    tcflush(STDIN_FILENO, TCIFLUSH);
#endif

    if (is_cancelled) {
        return 1;
    }

    std::cout << "\n=== LEONARD SELECTIONS ===\n";
    if (selected_sample) {
        std::cout << "Sample: " << selected_sample->name << "\n";
    }
    std::cout << "Variables: ";
    for (size_t i = 0; i < all_variables.size(); ++i) {
        if (var_states[i]) std::cout << all_variables[i] << " ";
    }
    std::cout << "\nPopulations: ";
    for (size_t i = 0; i < all_populations.size(); ++i) {
        if (pop_states[i]) std::cout << all_populations[i] << " ";
    }
    std::cout << "\nAnalysis Method:\n" << choices[analysis_choice] << "\n\n";
    return 0;
};
/*
Strip the .fcs from the file names
Rename Variables to Detectors, move populations over and insert a new one called Stains
not checkboxs just labels
When a user selects a sample fill in the stain column with $PnS from that sample
there is cruft before the first real population with a real checkbox.
It seems to be related to a label with special characters so they may need to be sanitized unless you have a better idea
probably have to be carful passing everything in quotes to command line
it leaves the screen in a nasty state lets do a screen clear on either escape or return.
Is it possible to change the bottom banner to include Space to select and Arrow keys to navigate?
I'm not sure the logic on the sample selection is right. Only populations that are in the first sample should be enabled.
Other opuplations that contain all the populations selected in the first sample can be selected but only if they contain all.
can we make tab break out of the selection menu boxes to the decision menu? 
*/