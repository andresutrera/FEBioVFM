// src/vfm/InternalWorkAssembler.hpp
#pragma once
#include <vector>
#include <string>
#include <functional>
#include <FECore/mat3d.h>
#include "build/mesh_info.hpp"
#include "domain/vfm_tensors.hpp"

class InternalWorkAssembler
{
public:
    using ParamSetter = std::function<bool(const std::vector<double> &, std::string &)>;
    using StressComputer = std::function<bool(std::string &)>; // fills Stresses (σ,P)
    using ToVirtGrad = std::function<mat3d(const mat3d &)>;    // VirtualGradientFromDeformation

    InternalWorkAssembler(const MeshDims &dims,
                          const MeshQuad &quad,
                          const VirtualDeformations &vdef,
                          const Stresses &stresses,
                          ParamSetter setParams,
                          StressComputer computeStress,
                          ToVirtGrad toVirtualGrad)
        : m_dims(dims), m_quad(quad), m_vdef(vdef), m_stress(stresses),
          m_setParams(std::move(setParams)),
          m_computeStress(std::move(computeStress)),
          m_toVG(std::move(toVirtualGrad)) {}

    // Returns flattened IW vector W[vf,t]
    std::vector<double> operator()(const std::vector<double> &params, std::string &err)
    {
        if (!m_setParams(params, err))
            return {};
        if (!m_computeStress(err))
            return {};

        const std::size_t VF = m_vdef.nVF();
        if (VF == 0)
            return {};
        // assume time alignment by index across P and virtual fields
        const std::size_t T = m_stress.nTimes();
        if (T == 0)
            return {};

        std::vector<double> W(VF * T, 0.0);

        for (std::size_t v = 0; v < VF; ++v)
        {
            const int vfTimes = m_vdef.nTimes((VFIdx)v);
            if (vfTimes == 0)
            {
                err = "virtual field has no time steps.";
                return {};
            }

            bool useSingleTime = false;
            if (vfTimes < static_cast<int>(T))
            {
                if (vfTimes == 1)
                {
                    useSingleTime = true;
                }
                else
                {
                    err = "virtual field has fewer time steps. A constant field can be defined using a single time entry.";
                    return {};
                }
            }
            for (std::size_t t = 0; t < T; ++t)
            {
                double acc = 0.0;

                // loop elements using precomputed gp layout
                const std::size_t nElem = m_quad.gpPerElem.size();
                for (std::size_t e = 0; e < nElem; ++e)
                {
                    const std::size_t off = m_quad.offset[e];
                    const std::size_t nint = m_quad.gpPerElem[e];

                    // optional: guard same shape
                    // if (m_stress.nGauss((TimeIdx)t, e) != nint) { err="shape mismatch"; return {}; }

                    for (std::size_t g = 0; g < nint; ++g)
                    {
                        const mat3d &P = m_stress.crefP((TimeIdx)t, e, g);
                        mat3d G;
                        if (useSingleTime)
                        {
                            G = m_toVG(m_vdef.crefF((VFIdx)v, (TimeIdx)0, e, g));
                        }
                        else
                        {
                            G = m_toVG(m_vdef.crefF((VFIdx)v, (TimeIdx)t, e, g));
                        }

                        acc += P.dotdot(G) * m_quad.jw[off + g]; // FEBio’s double contraction
                    }
                }
                W[v * T + t] = acc;
            }
        }
        return W;
    }

private:
    const MeshDims &m_dims;
    const MeshQuad &m_quad;
    const VirtualDeformations &m_vdef;
    const Stresses &m_stress;
    ParamSetter m_setParams;
    StressComputer m_computeStress;
    ToVirtGrad m_toVG;
};
