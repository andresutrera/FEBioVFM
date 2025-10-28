#pragma once
#include <string>
#include <vector>
#include <FECore/vec3d.h>



struct NodalSample { int id=-1; vec3d v; };                 // node id
struct TimeSliceNodes { int t=0; std::vector<NodalSample> nodes; };

struct SurfaceLoad { std::string surf; vec3d v; };
struct TimeSliceLoads { int t=0; std::vector<SurfaceLoad> loads; };

struct VirtualFieldXML {
    int id=-1;                               // VF id as integer; -1 if missing
    std::vector<TimeSliceNodes> times;       // per-time sparse nodal samples
};

struct XMLInput {
    std::vector<TimeSliceNodes>  measuredU;     // <MeasuredDisplacements>
    std::vector<VirtualFieldXML> virtualU;      // <VirtualDisplacements>
    std::vector<TimeSliceLoads>  measuredLoads; // <MeasuredLoads>
    struct Param { std::string name; double init=0, lo=0, hi=0, scale=1; };
    std::vector<Param> parameters;              // <Parameters>
    struct Options
    {
        enum class Type
        {
            Levmar,
            ConstrainedLevmar
        };

        struct Entry
        {
            bool set = false;
            double value = 0.0;
        };

        bool present = false;
        Type type = Type::ConstrainedLevmar;
        Entry tau;
        Entry gradTol;
        Entry stepTol;
        Entry objTol;
        Entry diffScale;
        Entry maxIters;
        bool planeDeformation = false;
        bool planeDeformationSet = false;
        std::string saveVirtualWork;
        bool saveVirtualWorkSet = false;
    } options;
};

class VFMXmlReader {
public:
    bool read(const char* path, XMLInput& out, std::string& err);
};
