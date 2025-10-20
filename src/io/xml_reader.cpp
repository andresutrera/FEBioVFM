#include "xml_reader.hpp"
#include <FEBioXML/xmltool.h>
#include <cstring>
#include <cstdlib>
#include <string>

static void ParseParameters(XMLTag& tag, XMLInput& out)
{
    ++tag;
    while (!tag.isend())
    {
        if (tag == "param")
        {
            XMLInput::Param p;
            if (const char* nm = tag.AttributeValue("name")) p.name = nm;
            double d[4] = {0,0,0,1};
            tag.value(d, 4);
            p.init  = d[0];
            p.lo    = d[1];
            p.hi    = d[2];
            p.scale = d[3];
            out.parameters.push_back(p);
        }
        else throw XMLReader::InvalidTag(tag);
        ++tag;
    }
}

// Parse <MeasuredDisplacements> or a single <virtualdisplacement> body
static void ParseDisplacementBlock(XMLTag& tag, std::vector<TimeSliceNodes>& dst)
{
    dst.clear();
    XMLTag timeTag(tag); ++timeTag;
    bool any = false;
    while (!timeTag.isend())
    {
        if (timeTag == "time")
        {
            any = true;

            // time attribute as integer
            int tval = 0;
            if (const char* ts = timeTag.AttributeValue("t", false)) tval = std::atoi(ts);

            TimeSliceNodes tsn; tsn.t = tval;

            XMLTag nodeTag = timeTag; ++nodeTag;
            while (!nodeTag.isend())
            {
                if ((nodeTag == "node") || (nodeTag == "elem"))
                {
                    const char* szId = nodeTag.AttributeValue("id");
                    int nid = szId ? std::atoi(szId) : -1;
                    double v[3] = {0,0,0};
                    nodeTag.value(v, 3);
                    tsn.nodes.push_back({nid, {v[0], v[1], v[2]}});
                }
                else throw XMLReader::InvalidTag(nodeTag);

                nodeTag.skip(); ++nodeTag;
            }

            dst.push_back(std::move(tsn));
        }
        timeTag.skip(); ++timeTag;
    }
    if (!any) throw XMLReader::InvalidTag(tag);
}

static void ParseMeasuredLoads(XMLTag& tag, std::vector<TimeSliceLoads>& dst)
{
    dst.clear();
    XMLTag timeTag(tag); ++timeTag;
    bool any = false;
    while (!timeTag.isend())
    {
        if (timeTag == "time")
        {
            any = true;

            int tval = 0;
            if (const char* ts = timeTag.AttributeValue("t", false)) tval = std::atoi(ts);

            TimeSliceLoads tsl; tsl.t = tval;

            XMLTag surfTag = timeTag; ++surfTag;
            while (!surfTag.isend())
            {
                if (surfTag == "surface")
                {
                    const char* sid = surfTag.AttributeValue("id");
                    std::string s = sid ? sid : "";
                    double f[3] = {0,0,0};
                    surfTag.value(f, 3);
                    tsl.loads.push_back({s, {f[0], f[1], f[2]}});
                }
                else throw XMLReader::InvalidTag(surfTag);

                surfTag.skip(); ++surfTag;
            }

            dst.push_back(std::move(tsl));
        }
        timeTag.skip(); ++timeTag;
    }
    if (!any) throw XMLReader::InvalidTag(tag);
}

bool VFMXmlReader::read(const char* path, XMLInput& out, std::string& err)
{
    out = XMLInput{};

    XMLReader xml;
    if (!xml.Open(path)) { err = "Failed to open XML file."; return false; }

    XMLTag root;
    if (!xml.FindTag("febio_optimize", root)) { err = "Missing <febio_optimize> root."; xml.Close(); return false; }

    const char* ver = root.AttributeValue("version", true);
    if ((ver == nullptr) || (std::strcmp(ver, "2.0") != 0)) {
        err = "Invalid <febio_optimize> version. Expected 2.0.";
        xml.Close(); return false;
    }

    bool ok = true;
    try {
        ++root;
        while (!root.isend())
        {
            if (root.Name()[0] == '#') { root.skip(); ++root; continue; }

            if (root == "Parameters") {
                ParseParameters(root, out);
                root.skip();
            }
            else if (root == "MeasuredDisplacements") {
                ParseDisplacementBlock(root, out.measuredU);
                root.skip();
            }
            else if (root == "VirtualDisplacements") {
                bool foundAny = false;

                XMLTag vfTag(root); ++vfTag;
                while (!vfTag.isend())
                {
                    if (vfTag == "virtualdisplacement")
                    {
                        foundAny = true;
                        VirtualFieldXML vf;

                        // VF id as integer; default -1 if missing
                        if (const char* id = vfTag.AttributeValue("id", false))
                            vf.id = std::atoi(id);

                        ParseDisplacementBlock(vfTag, vf.times);
                        out.virtualU.push_back(std::move(vf));
                    }
                    else if (vfTag == "time")
                    {
                        // Legacy mode: treat as one anonymous VF with id = -1
                        foundAny = true;
                        VirtualFieldXML vf; vf.id = -1;
                        ParseDisplacementBlock(root, vf.times);
                        out.virtualU.push_back(std::move(vf));
                        break;
                    }
                    else throw XMLReader::InvalidTag(vfTag);

                    vfTag.skip(); ++vfTag;
                }
                if (!foundAny) throw XMLReader::InvalidTag(root);
                root.skip();
            }
            else if (root == "MeasuredLoads") {
                ParseMeasuredLoads(root, out.measuredLoads);
                root.skip();
            }
            else {
                throw XMLReader::InvalidTag(root);
            }

            ++root;
        }
    }
    catch (const XMLReader::Error&) {
        err = "XML parse error in VFM input.";
        ok = false;
    }
    catch (...) {
        err = "Unknown exception in VFM XML reader.";
        ok = false;
    }

    xml.Close();
    return ok;
}
