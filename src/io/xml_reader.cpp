#include "xml_reader.hpp"
#include <FEBioXML/xmltool.h>
#include <cstring>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <cctype>

static bool parse_bool_tag_value(XMLTag &tag)
{
    const char *sz = tag.szvalue();
    if (sz == nullptr)
        throw XMLReader::InvalidTag(tag);

    std::string text(sz);
    auto notSpace = [](unsigned char ch)
    { return std::isspace(ch) == 0; };

    auto first = std::find_if(text.begin(), text.end(), notSpace);
    auto last = std::find_if(text.rbegin(), text.rend(), notSpace).base();
    if (first < last)
        text.assign(first, last);
    else
        text.clear();

    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c)
                   { return static_cast<char>(std::tolower(c)); });

    if (text == "true" || text == "1")
        return true;
    if (text == "false" || text == "0")
        return false;

    throw XMLReader::InvalidTag(tag);
}

static std::string parse_trimmed_tag_value(XMLTag &tag)
{
    const char *sz = tag.szvalue();
    if (sz == nullptr)
        throw XMLReader::InvalidTag(tag);

    std::string text(sz);
    auto notSpace = [](unsigned char ch)
    { return std::isspace(ch) == 0; };

    auto first = std::find_if(text.begin(), text.end(), notSpace);
    auto last = std::find_if(text.rbegin(), text.rend(), notSpace).base();
    if (first < last)
        text.assign(first, last);
    else
        text.clear();

    return text;
}

static void ParseParameters(XMLTag &tag, XMLInput &out)
{
    ++tag;
    while (!tag.isend())
    {
        if (tag == "param")
        {
            XMLInput::Param p;
            if (const char *nm = tag.AttributeValue("name"))
                p.name = nm;
            double d[4] = {0, 0, 0, 1};
            tag.value(d, 4);
            p.init = d[0];
            p.lo = d[1];
            p.hi = d[2];
            p.scale = d[3];
            out.parameters.push_back(p);
        }
        else
            throw XMLReader::InvalidTag(tag);
        ++tag;
    }
}

static void ParseGeneralOptions(XMLTag &tag, XMLInput &out)
{
    XMLTag child(tag);
    ++child;
    while (!child.isend())
    {
        if (child == "plane_deformation")
        {
            const bool value = parse_bool_tag_value(child);
            out.options.planeDeformation = value;
            out.options.planeDeformationSet = true;
        }
        else if (child == "save_virtual_work")
        {
            const std::string value = parse_trimmed_tag_value(child);
            std::string lower = value;
            std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c)
                           { return static_cast<char>(std::tolower(c)); });
            if (!lower.empty() && lower == "false")
            {
                out.options.saveVirtualWork.clear();
                out.options.saveVirtualWorkSet = false;
            }
            else
            {
                out.options.saveVirtualWork = value;
                out.options.saveVirtualWorkSet = !value.empty();
            }
        }
        else
            throw XMLReader::InvalidTag(child);

        child.skip();
        ++child;
    }
}

static void ParseOptions(XMLTag &tag, XMLInput &out)
{
    XMLInput::Options opts = out.options;
    opts.present = true;

    if (const char *typeAttr = tag.AttributeValue("type", false))
    {
        std::string type = typeAttr;
        std::transform(type.begin(), type.end(), type.begin(), [](unsigned char c)
                       { return static_cast<char>(std::tolower(c)); });
        if (type == "levmar")
            opts.type = XMLInput::Options::Type::Levmar;
        else if (type == "constrained levmar")
            opts.type = XMLInput::Options::Type::ConstrainedLevmar;
        else
            throw XMLReader::InvalidTag(tag);
    }

    XMLTag child(tag);
    ++child;
    while (!child.isend())
    {
        double val[1] = {0.0};
        if (child == "tau")
        {
            child.value(val, 1);
            opts.tau.set = true;
            opts.tau.value = val[0];
        }
        else if (child == "grad_tol")
        {
            child.value(val, 1);
            opts.gradTol.set = true;
            opts.gradTol.value = val[0];
        }
        else if (child == "step_tol")
        {
            child.value(val, 1);
            opts.stepTol.set = true;
            opts.stepTol.value = val[0];
        }
        else if (child == "obj_tol")
        {
            child.value(val, 1);
            opts.objTol.set = true;
            opts.objTol.value = val[0];
        }
        else if (child == "f_diff_scale")
        {
            child.value(val, 1);
            opts.diffScale.set = true;
            opts.diffScale.value = val[0];
        }
        else if (child == "max_iter")
        {
            child.value(val, 1);
            opts.maxIters.set = true;
            opts.maxIters.value = val[0];
        }
        else
            throw XMLReader::InvalidTag(child);

        child.skip();
        ++child;
    }

    out.options = opts;
}

// Parse <MeasuredDisplacements> or a single <virtualdisplacement> body
static void ParseDisplacementBlock(XMLTag &tag, std::vector<TimeSliceNodes> &dst)
{
    dst.clear();
    XMLTag timeTag(tag);
    ++timeTag;
    bool any = false;
    while (!timeTag.isend())
    {
        if (timeTag == "time")
        {
            any = true;

            // time attribute as integer
            int tval = 0;
            if (const char *ts = timeTag.AttributeValue("t", false))
                tval = std::atoi(ts);

            TimeSliceNodes tsn;
            tsn.t = tval;

            XMLTag nodeTag = timeTag;
            ++nodeTag;
            while (!nodeTag.isend())
            {
                if ((nodeTag == "node") || (nodeTag == "elem"))
                {
                    const char *szId = nodeTag.AttributeValue("id");
                    int nid = szId ? std::atoi(szId) : -1;
                    double v[3] = {0, 0, 0};
                    nodeTag.value(v, 3);
                    tsn.nodes.push_back({nid, {v[0], v[1], v[2]}});
                }
                else
                    throw XMLReader::InvalidTag(nodeTag);

                nodeTag.skip();
                ++nodeTag;
            }

            dst.push_back(std::move(tsn));
        }
        timeTag.skip();
        ++timeTag;
    }
    if (!any)
        throw XMLReader::InvalidTag(tag);
}

static void ParseMeasuredLoads(XMLTag &tag, std::vector<TimeSliceLoads> &dst)
{
    dst.clear();
    XMLTag timeTag(tag);
    ++timeTag;
    bool any = false;
    while (!timeTag.isend())
    {
        if (timeTag == "time")
        {
            any = true;

            int tval = 0;
            if (const char *ts = timeTag.AttributeValue("t", false))
                tval = std::atoi(ts);

            TimeSliceLoads tsl;
            tsl.t = tval;

            XMLTag surfTag = timeTag;
            ++surfTag;
            while (!surfTag.isend())
            {
                if (surfTag == "surface")
                {
                    const char *sid = surfTag.AttributeValue("id");
                    std::string s = sid ? sid : "";
                    double f[3] = {0, 0, 0};
                    surfTag.value(f, 3);
                    tsl.loads.push_back({s, {f[0], f[1], f[2]}});
                }
                else
                    throw XMLReader::InvalidTag(surfTag);

                surfTag.skip();
                ++surfTag;
            }

            dst.push_back(std::move(tsl));
        }
        timeTag.skip();
        ++timeTag;
    }
    if (!any)
        throw XMLReader::InvalidTag(tag);
}

bool VFMXmlReader::read(const char *path, XMLInput &out, std::string &err)
{
    out = XMLInput{};

    XMLReader xml;
    if (!xml.Open(path))
    {
        err = "Failed to open XML file.";
        return false;
    }

    XMLTag root;
    if (!xml.FindTag("febio_optimize", root))
    {
        err = "Missing <febio_optimize> root.";
        xml.Close();
        return false;
    }

    const char *ver = root.AttributeValue("version", true);
    if ((ver == nullptr) || (std::strcmp(ver, "2.0") != 0))
    {
        err = "Invalid <febio_optimize> version. Expected 2.0.";
        xml.Close();
        return false;
    }

    bool ok = true;
    try
    {
        ++root;
        while (!root.isend())
        {
            if (root.Name()[0] == '#')
            {
                root.skip();
                ++root;
                continue;
            }

            if (root == "Parameters")
            {
                ParseParameters(root, out);
                root.skip();
            }
            else if (root == "Options")
            {
                ParseGeneralOptions(root, out);
                root.skip();
            }
            else if (root == "MeasuredDisplacements")
            {
                ParseDisplacementBlock(root, out.measuredU);
                root.skip();
            }
            else if (root == "VirtualDisplacements")
            {
                bool foundAny = false;

                XMLTag vfTag(root);
                ++vfTag;
                while (!vfTag.isend())
                {
                    if (vfTag == "virtualdisplacement")
                    {
                        foundAny = true;
                        VirtualFieldXML vf;

                        // VF id as integer; default -1 if missing
                        if (const char *id = vfTag.AttributeValue("id", false))
                            vf.id = std::atoi(id);

                        ParseDisplacementBlock(vfTag, vf.times);
                        out.virtualU.push_back(std::move(vf));
                    }
                    else if (vfTag == "time")
                    {
                        // Legacy mode: treat as one anonymous VF with id = -1
                        foundAny = true;
                        VirtualFieldXML vf;
                        vf.id = -1;
                        ParseDisplacementBlock(root, vf.times);
                        out.virtualU.push_back(std::move(vf));
                        break;
                    }
                    else
                        throw XMLReader::InvalidTag(vfTag);

                    vfTag.skip();
                    ++vfTag;
                }
                if (!foundAny)
                    throw XMLReader::InvalidTag(root);
                root.skip();
            }
            else if (root == "MeasuredLoads")
            {
                ParseMeasuredLoads(root, out.measuredLoads);
                root.skip();
            }
            else if (root == "Optimization")
            {
                ParseOptions(root, out);
                root.skip();
            }
            else
            {
                throw XMLReader::InvalidTag(root);
            }

            ++root;
        }
    }
    catch (const XMLReader::Error &)
    {
        err = "XML parse error in VFM input.";
        ok = false;
    }
    catch (...)
    {
        err = "Unknown exception in VFM XML reader.";
        ok = false;
    }

    xml.Close();
    return ok;
}
