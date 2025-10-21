#include "optimization/vfm_problem.hpp"
#include <algorithm>
#include <vector>
#include <FECore/FEModel.h>
#include <FECore/FEMaterial.h>
#include <FECore/log.h>
#include "io/loader.hpp"
#include "FE/params_febio.hpp"
#include "FE/shape_provider_febio.hpp"
#include "services/kinematics.hpp"
#include "FE/material_provider_febio.hpp"
#include "services/stress_eval.hpp"
#include "diag/felog_bridge.hpp"

void VFMProblem::reset()
{
  fem = nullptr;
  dims = MeshDims{};
  conn = MeshConn{};
  quad = MeshQuad{};
  surfaces.clear();
  state.clear();
}

namespace
{
  void collect_surface_names(const MeasuredLoad &loads, std::vector<std::string> &out)
  {
    out.clear();
    const std::size_t T = loads.nTimes();
    for (std::size_t t = 0; t < T; ++t)
    {
      const LoadFrame &frame = loads.frame(static_cast<TimeIdx>(t));
      for (const auto &entry : frame.loads)
      {
        if (entry.surface.empty())
          continue;
        if (std::find(out.begin(), out.end(), entry.surface) == out.end())
          out.push_back(entry.surface);
      }
    }
  }
} // namespace

bool prepare_vfm_problem(FEModel &fem,
                         const XMLInput &input,
                         VFMProblem &problem,
                         std::string &err)
{
  diag::ScopedFEBind bind(&fem);
  problem.reset();
  problem.fem = &fem;

  // mesh info
  if (!build_mesh_info(fem, problem.dims, problem.conn, problem.quad, err))
    return false;

  // load data from XML DTO
  MeasuredData measured;
  VirtualFields virtuals;
  MeasuredLoad loads;

  if (!VFMLoader::load_measuredU(input, problem.dims, measured, err))
    return false;
  if (!VFMLoader::load_virtualU(input, problem.dims, virtuals, err))
    return false;
  if (!VFMLoader::load_measuredF(input, problem.dims, loads, err))
    return false;
  if (!VFMLoader::load_params(input, problem.state, err))
    return false;

  feLog("Success loading input data.\n");

  // move data into state and configure tensors
  problem.state.measured = std::move(measured);
  problem.state.virtuals = std::move(virtuals);
  problem.state.loads = std::move(loads);
  problem.state.configure_tensors(problem.quad.gpPerElem, problem.state.virtuals.nVF());
  problem.state.mirror_frames_from_displacements();

  // kinematics
  FEBioShapeProvider shp(problem.conn);
  std::string kinErr;
  Kinematics::compute_measured(problem.quad, shp, problem.state.measured, problem.state.def, true, kinErr);
  Kinematics::compute_virtuals(problem.quad, shp, problem.state.virtuals, problem.state.vdef, true, kinErr);
  if (!kinErr.empty())
  {
    err = kinErr;
    return false;
  }

  feLog("Success computed VFM kinetics.\n");

  // precompute surface info
  std::vector<std::string> surfaceNames;
  collect_surface_names(problem.state.loads, surfaceNames);
  if (!surfaceNames.empty())
  {
    if (!build_surface_info(fem.GetMesh(), problem.dims, surfaceNames, problem.surfaces, err))
      return false;
  }

  return true;
}
