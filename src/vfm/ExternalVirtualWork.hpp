#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <FECore/FEModel.h>
#include <FECore/vec3d.h>
#include "build/mesh_info.hpp"
#include "build/surface_info.hpp"
#include "domain/vfm_displacements.hpp"

class ExternalVirtualWorkAssembler
{
public:
    ExternalVirtualWorkAssembler(FEModel &fem,
                                 const MeshDims &dims,
                                 const VirtualFields &virtuals,
                                 const MeasuredLoad &loads)
        : m_fem(fem), m_dims(dims), m_virtuals(virtuals), m_loads(loads) {}

    // Returns flattened EW vector W[vf, t]
    std::vector<double> operator()(std::string &err)
    {
        err.clear();

        const std::size_t VF = m_virtuals.nVF();
        const std::size_t T = m_loads.nTimes();
        if (VF == 0 || T == 0)
            return {};

        std::vector<std::string> surfaceNames;
        for (std::size_t t = 0; t < T; ++t)
        {
            const LoadFrame &frame = m_loads.frame(static_cast<TimeIdx>(t));
            for (const auto &entry : frame.loads)
            {
                if (std::find(surfaceNames.begin(), surfaceNames.end(), entry.surface) == surfaceNames.end())
                    surfaceNames.push_back(entry.surface);
            }
        }

        SurfaceMap surfaces;
        if (!surfaceNames.empty())
        {
            if (!build_surface_info(m_fem.GetMesh(), m_dims, surfaceNames, surfaces, err))
                return {};
        }

        std::vector<double> W(VF * T, 0.0);

        for (std::size_t v = 0; v < VF; ++v)
        {
            const auto &vfSeries = m_virtuals.getVF(static_cast<VFIdx>(v));
            if (vfSeries.nTimes() < T)
            {
                err = "virtual field has fewer time steps than loads";
                return {};
            }

            for (std::size_t t = 0; t < T; ++t)
            {
                const LoadFrame &frame = m_loads.frame(static_cast<TimeIdx>(t));
                const VirtualFrame &vfFrame = vfSeries.getTime(static_cast<TimeIdx>(t));
                double acc = 0.0;

                for (const auto &entry : frame.loads)
                {
                    auto it = surfaces.find(entry.surface);
                    if (it == surfaces.end())
                    {
                        err = "missing surface mapping for " + entry.surface;
                        return {};
                    }

                    const auto &nodes = it->second.idx;
                    if (nodes.empty())
                    {
                        err = "surface with no nodes: " + entry.surface;
                        return {};
                    }

                    vec3d ustar(0.0, 0.0, 0.0);
                    bool hasNode = false;
                    for (std::size_t idx : nodes)
                    {
                        ustar = vfFrame.u.getNode(idx);
                        hasNode = true;
                        break;
                    }
                    if (!hasNode)
                    {
                        err = "surface nodes unavailable for " + entry.surface;
                        return {};
                    }

                    acc += entry.force * ustar;
                }

                W[v * T + t] = acc;
            }
        }

        return W;
    }

private:
    FEModel &m_fem;
    const MeshDims &m_dims;
    const VirtualFields &m_virtuals;
    const MeasuredLoad &m_loads;
};
