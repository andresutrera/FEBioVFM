#include "vfm_task.h"
#include "io/loader.hpp"
#include "domain/vfm_displacements.hpp"
#include "state/vfm_state.hpp"

#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FEMaterial.h>

#include "FE/shape_provider_febio.hpp"
#include "services/kinematics.hpp"

#include "FE/material_provider_febio.hpp"
#include "services/stress_eval.hpp"

#include "FE/params_febio.hpp"
#include "vfm/InternalWork.hpp"
#include "vfm/ExternalVirtualWork.hpp"
#include "optimization/vfm_solver.hpp"
#include <cmath>

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

    std::string err;
    VFMXmlReader rdr;
    if (!rdr.read(m_inputPath.c_str(), m_input, err))
    {
        feLogError(err.c_str());
        return false;
    }

    m_state.clear();

    MeshDims dims;
    MeshConn conn;
    MeshQuad quad;
    if (!build_mesh_info(*GetFEModel(), dims, conn, quad, err))
    {
        feLogError(err.c_str());
        return false;
    }

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

    m_dims = std::move(dims);
    m_conn = std::move(conn);
    m_quad = std::move(quad);

    // move into state BEFORE configuring tensors
    m_state.measured = std::move(measured);
    m_state.virtuals = std::move(virtuals);
    m_state.loads = std::move(loads);
    // 4) Configure tensor containers (shape + nVF)
    m_state.configure_tensors(m_quad.gpPerElem, m_state.virtuals.nVF());

    // 5) Create tensor frames mirroring displacement timelines
    m_state.mirror_frames_from_displacements();

    FEBioShapeProvider shp(m_conn);
    Kinematics::compute_measured(m_quad, shp, m_state.measured, m_state.def, true, err);
    Kinematics::compute_virtuals(m_quad, shp, m_state.virtuals, m_state.vdef, true, err);

    // after compute_*
    if (!err.empty())
    {
        feLogError(err.c_str());
        return false;
    }

    FEModel &fem = *GetFEModel();
    FEBioParameterApplier setParams(fem, m_state);

    // force a change
    std::vector<double> p(m_state.params.size());
    for (size_t i = 0; i < p.size(); ++i)
        p[i] = m_state.params[i].value;

    if (!setParams.apply(p, err))
    {
        feLogError(err.c_str());
        return false;
    }

    // ensure fresh model state
    for (int i = 0; i < fem.Materials(); ++i)
        fem.GetMaterial(i)->Init();

    FEBioMaterialProvider mat(m_conn);

    if (!StressEval::cauchy(m_state.def, m_state.stresses, mat, err))
    {
        feLogError(err.c_str());
        return false;
    }

    if (!StressEval::first_piola(m_state.def, m_state.stresses, m_state.stresses, err))
    {
        feLogError(err.c_str());
        return false;
    }

    return true;
}

bool VFMTask::Run()
{
    feLog("...........................................................................\n");
    feLog("                                    RUN                                    \n");
    feLog("...........................................................................\n");
    feLog("\n");
    feLog("\n");
    std::string err;

    if (m_state.params.empty())
    {
        feLog("No parameters to optimize.\n");
        return true;
    }

    FEModel &fem = *GetFEModel();
    FEBioParameterApplier setParams(fem, m_state);
    FEBioMaterialProvider mat(m_conn);

    auto paramSetter = [&](const std::vector<double> &values, std::string &perr) -> bool
    {
        return setParams.apply(values, perr);
    };

    auto computeStress = [&](std::string &perr) -> bool
    {
        if (!StressEval::cauchy(m_state.def, m_state.stresses, mat, perr))
            return false;
        return StressEval::first_piola(m_state.def, m_state.stresses, m_state.stresses, perr);
    };

    auto toVirtGrad = [](const mat3d &Fstar)
    {
        mat3d G = Fstar;
        G[0][0] -= 1.0;
        G[1][1] -= 1.0;
        G[2][2] -= 1.0;
        return G;
    };

    InternalWorkAssembler internal(m_dims, m_quad, m_state.vdef, m_state.stresses,
                                   paramSetter, computeStress, toVirtGrad);

    ExternalVirtualWorkAssembler external(fem, m_dims, m_state.virtuals, m_state.loads);
    std::vector<double> ew = external(err);
    if (!err.empty())
    {
        feLogError(err.c_str());
        return false;
    }
    if (ew.empty())
    {
        feLog("External work vector empty. Nothing to optimize.\n");
        return true;
    }

    std::vector<double> params;
    params.reserve(m_state.params.size());
    for (const auto &q : m_state.params)
        params.push_back(q.value);

    if (!run_vfm_levmar(params, internal, ew, 100, err))
    {
        if (!err.empty())
            feLogError(err.c_str());
        else
            feLogError("Levenberg–Marquardt solver failed.");
        return false;
    }

    if (!setParams.apply(params, err))
    {
        feLogError(err.c_str());
        return false;
    }

    feLog("Optimized parameters:\n");
    for (size_t i = 0; i < m_state.params.size(); ++i)
    {
        m_state.params[i].value = params[i];
        feLog("  %s = %.6g\n", m_state.params[i].spec.name.c_str(), params[i]);
    }

    // Refresh stress state with optimized parameters for downstream use.
    if (!computeStress(err))
    {
        feLogError(err.c_str());
        return false;
    }

    feLog("Optimization complete.\n");
    return true;
}
