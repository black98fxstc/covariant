#include <iostream>
#include <string>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <libxml/tree.h>

void inject_gating_ml(const std::string& wsp_file, const std::string& gml_file, const std::string& out_wsp_file) {
    xmlInitParser();

    xmlDocPtr wspDoc = xmlParseFile(wsp_file.c_str());
    xmlDocPtr gmlDoc = xmlParseFile(gml_file.c_str());

    if (!wspDoc || !gmlDoc) {
        std::cerr << "Error: Could not parse input files.\n";
        return;
    }

    // Locate the first Subpopulations node under a SampleNode to inject into.
    // (In a more advanced implementation, you could target a specific Sample by name/ID)
    xmlXPathContextPtr wspCtx = xmlXPathNewContext(wspDoc);
    xmlXPathObjectPtr sampleNodeObj = xmlXPathEvalExpression((const xmlChar*)"(//SampleNode/Subpopulations)[1]", wspCtx);
    
    if (sampleNodeObj && sampleNodeObj->nodesetval && sampleNodeObj->nodesetval->nodeNr > 0) {
        xmlNodePtr subpopsNode = sampleNodeObj->nodesetval->nodeTab[0];

        // Find all Gates inside the incoming Gating-ML document
        xmlXPathContextPtr gmlCtx = xmlXPathNewContext(gmlDoc);
        xmlXPathObjectPtr gatesObj = xmlXPathEvalExpression((const xmlChar*)"//*[contains(local-name(), 'Gate')]", gmlCtx);

        if (gatesObj && gatesObj->nodesetval) {
            for (int i = 0; i < gatesObj->nodesetval->nodeNr; ++i) {
                xmlNodePtr gateNode = gatesObj->nodesetval->nodeTab[i];
                
                // Wrap the Gating-ML element in FlowJo's expected <Population><Gate>...</Gate></Population>
                xmlNodePtr popNode = xmlNewNode(NULL, (const xmlChar*)"Population");
                std::string popName = "Imported_Alg_Gate_" + std::to_string(i + 1);
                
                xmlSetProp(popNode, (const xmlChar*)"name", (const xmlChar*)popName.c_str());
                xmlSetProp(popNode, (const xmlChar*)"expanded", (const xmlChar*)"1");
                xmlSetProp(popNode, (const xmlChar*)"owningGroup", (const xmlChar*)"");
                
                xmlNodePtr wrapperGateNode = xmlNewNode(NULL, (const xmlChar*)"Gate");
                xmlAddChild(popNode, wrapperGateNode);
                
                // Deep copy the incoming Gating-ML gate node
                xmlNodePtr copiedGate = xmlCopyNode(gateNode, 1);
                xmlAddChild(wrapperGateNode, copiedGate);
                
                // Append to the Workspace Subpopulations
                xmlAddChild(subpopsNode, popNode);
            }
        }
        xmlXPathFreeObject(gatesObj);
        xmlXPathFreeContext(gmlCtx);
        
        xmlSaveFormatFileEnc(out_wsp_file.c_str(), wspDoc, "UTF-8", 1);
        std::cout << "Successfully injected Gating-ML and saved to: " << out_wsp_file << "\n";
    } else {
        std::cerr << "Failed to find a Subpopulations node in the workspace to inject into.\n";
    }

    xmlXPathFreeObject(sampleNodeObj);
    xmlXPathFreeContext(wspCtx);
    xmlFreeDoc(gmlDoc);
    xmlFreeDoc(wspDoc);
    xmlCleanupParser();
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <in.wsp> <in.gml> <out.wsp>\n";
        return 1;
    }

    inject_gating_ml(argv[1], argv[2], argv[3]);
    return 0;
}