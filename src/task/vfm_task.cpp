#include "vfm_task.h"
#include "build/mesh_info.hpp"
#include "io/loader.hpp"
#include "domain/vfm_displacements.hpp"
#include "state/vfm_state.hpp"

#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>

#include "FE/shape_provider_febio.hpp"
#include "services/kinematics.hpp"

#include "FE/material_provider_febio.hpp"
#include "services/stress_eval.hpp"

VFMTask::VFMTask(FEModel *fem) : FECoreTask(fem) {}

bool VFMTask::Init(const char *xmlPath)
{
    m_inputPath = (xmlPath && *xmlPath) ? xmlPath : std::string{};
    feLog("\n");
    feLog("===========================================================================\n");
    feLog("                        VIRTUAL FIELDS METHOD (VFM)                        \n");
    feLog("===========================================================================\n");
    feLog("\n");
    feLog("...........................................................................\n");
    feLog("                                   SETUP                                   \n");
    feLog("...........................................................................\n");
    feLog("\n");
    feLog("\n");

    // parse XML → DTO
    std::string err;
    VFMXmlReader rdr;
    if (!rdr.read(m_inputPath.c_str(), m_input, err))
    {
        feLogError(err.c_str());
        return false;
    }

    // minimal report
    feLog("XML parsed: %s\n", m_inputPath.c_str());
    feLog("  measuredU: %zu time slices\n", m_input.measuredU.size());
    feLog("  virtualU : %zu fields\n", m_input.virtualU.size());
    feLog("  loads    : %zu time slices\n", m_input.measuredLoads.size());
    feLog("  params   : %zu\n", m_input.parameters.size());

    MeshDims dims;
    MeshConn conn;
    MeshQuad quad;
    if (!build_mesh_info(*GetFEModel(), dims, conn, quad, err))
    {
        feLogError(err.c_str());
        return false;
    }
    feLog("mesh: nodes=%zu elems=%zu\n", dims.nNodes, dims.nElems);

    MeasuredData measured;
    VirtualFields virtuals;
    MeasuredLoad loads;

    if (!VFMLoader::load_measuredU(m_input, dims, measured, err))
    {
        feLogError(err.c_str());
        return false;
    }
    if (!VFMLoader::load_virtualU(m_input, dims, virtuals, err))
    {
        feLogError(err.c_str());
        return false;
    }
    if (!VFMLoader::load_measuredF(m_input, dims, loads, err))
    {
        feLogError(err.c_str());
        return false;
    }
    if (!VFMLoader::load_params(m_input, m_state, err))
    {
        feLogError(err.c_str());
        return false;
    }

    feLog("VFM: params: %zu\n", m_state.params.size());
    for (const auto &q : m_state.params)
        feLog("  %s = %.6g  [%g, %g]  scale=%g\n",
              q.spec.name.c_str(), q.value, q.spec.lo, q.spec.hi, q.spec.scale);

    // move into state BEFORE configuring tensors
    m_state.clear();
    m_state.measured = std::move(measured);
    m_state.virtuals = std::move(virtuals);
    m_state.loads = std::move(loads);
    // 4) Configure tensor containers (shape + nVF)
    m_state.configure_tensors(quad.gpPerElem, m_state.virtuals.nVF());

    // 5) Create tensor frames mirroring displacement timelines
    m_state.mirror_frames_from_displacements();

    // Quick sizes
    feLog("mesh nodes=%zu elems=%zu gp[0]=%zu\n",
          dims.nNodes, dims.nElems, quad.gpPerElem.empty() ? 0 : quad.gpPerElem[0]);
    feLog("measuredU times=%zu  virtualU nVF=%zu\n",
          m_state.measured.series.nTimes(), m_state.virtuals.nVF());

    FEBioShapeProvider shp(conn);
    Kinematics::compute_measured(quad, shp, m_state.measured, m_state.def, true, err);
    Kinematics::compute_virtuals(quad, shp, m_state.virtuals, m_state.vdef, true, err);

    // after compute_*
    if (!err.empty())
    {
        feLogError(err.c_str());
        return false;
    }

    // measured F sizes
    feLog("F(meas): frames=%zu\n", m_state.def.nTimes());
    if (m_state.def.nTimes() > 0)
    {
        feLog("  elems=%zu  gp(e0)=%zu\n",
              m_state.def.nElements(0),
              m_state.def.nGauss(0, 0));
        // sample detF
        const mat3d &F00 = m_state.def.crefF(0, 0, 0);
        feLog("  detF[t0,e0,g0]=%g\n", F00.det());
    }

    // virtual F sizes
    feLog("F(virt): nVF=%zu\n", m_state.vdef.nVF());
    if (m_state.vdef.nVF() > 0)
    {
        feLog("  vf0.frames=%zu\n", m_state.vdef.nTimes(0));
        if (m_state.vdef.nTimes(0) > 0)
        {
            feLog("  vf0.elems=%zu  gp(e0)=%zu\n",
                  m_state.vdef.nElements(0, 0),
                  m_state.vdef.nGauss(0, 0, 0));
            const mat3d &Fv000 = m_state.vdef.crefF(0, 0, 0, 0);
            feLog("  detFv[v0,t0,e0,g0]=%g\n", Fv000.det());
        }
    }

    FEBioMaterialProvider mat(conn);

    m_state.stresses.addTime(); // one frame per def frame
    if (!StressEval::cauchy(m_state.def, m_state.stresses, mat, err))
    {
        feLogError(err.c_str());
        return false;
    }

    Stresses P;
    P.setElemShape(quad.gpPerElem);
    for (size_t k = 0; k < m_state.def.nTimes(); ++k)
        (void)P.addTime();
    if (!StressEval::first_piola(m_state.def, m_state.stresses, P, err))
    {
        feLogError(err.c_str());
        return false;
    }
    // store P if you keep it, or extend Stresses to hold both σ and P as you already do.

    feLog("stress: frames=%zu\n", m_state.stresses.nTimes());
    if (m_state.stresses.nTimes() > 0)
    {
        feLog("  elems=%zu  gp(e0)=%zu\n",
              m_state.stresses.nElements(0),
              m_state.stresses.nGauss(0, 0));
        const mat3d &s00 = m_state.stresses.crefSigma(10, 0, 0);
        const mat3d &p00 = P.crefP(10, 0, 0); // from the P container above
        feLog("  sigma00=%g  P00=%g\n", s00[2][2], p00[0][0]);
        feLog("  sigma [%g %g %g]\n", s00[0][0], s00[0][1], s00[0][2]);
        feLog("  sigma [%g %g %g]\n", s00[1][0], s00[1][1], s00[1][2]);
        feLog("  sigma [%g %g %g]\n", s00[2][0], s00[2][1], s00[2][2]);
    }

    return true;
}

bool VFMTask::Run()
{
    feLog("\n=== VFM RUN ===\n");
    feLog("Nothing to do yet.\n");
    return true;
}
