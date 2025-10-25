#pragma once
#include <vector>
#include <string>
#include <FECore/vec3d.h>
#include "build/surface_info.hpp"
#include "domain/vfm_displacements.hpp"

class ExternalVirtualWorkAssembler
{
public:
    ExternalVirtualWorkAssembler(const SurfaceMap &surfaces,
                                 const VirtualFields &virtuals,
                                 const MeasuredLoad &loads)
        : m_surfaces(surfaces), m_virtuals(virtuals), m_loads(loads) {}

    // Returns flattened EW vector W[vf, t]
    std::vector<double> operator()(std::string &err)
    {
        err.clear();

        const std::size_t VF = m_virtuals.nVF();
        const std::size_t T = m_loads.nTimes();

        bool useSingleTime = false;
        if (VF == 0 || T == 0)
            return {};

        std::vector<double> W(VF * T, 0.0);

        for (std::size_t v = 0; v < VF; ++v)
        {
            const auto &vfSeries = m_virtuals.getVF(static_cast<VFIdx>(v));
            if (vfSeries.nTimes() < T)
            {
                // err = "virtual field has fewer time steps than loads";
                // return {};
                useSingleTime = true;
            }

            for (std::size_t t = 0; t < T; ++t)
            {
                VirtualFrame vfFrame;
                const LoadFrame &frame = m_loads.frame(static_cast<TimeIdx>(t));
                if (useSingleTime)
                {
                    vfFrame = vfSeries.getTime(static_cast<TimeIdx>(0));
                }
                else
                {
                    vfFrame = vfSeries.getTime(static_cast<TimeIdx>(t));
                }

                double acc = 0.0;

                for (const auto &entry : frame.loads)
                {
                    auto it = m_surfaces.find(entry.surface);
                    if (it == m_surfaces.end())
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
    const SurfaceMap &m_surfaces;
    const VirtualFields &m_virtuals;
    const MeasuredLoad &m_loads;
};
