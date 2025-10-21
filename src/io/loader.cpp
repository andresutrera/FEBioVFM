#include "io/loader.hpp"
#include "state/vfm_state.hpp"
#include <cstdio>

namespace
{
  inline bool map_node(const MeshDims &dims, int nodeId, std::size_t &idx)
  {
    auto it = dims.nodeId2idx.find(nodeId);
    if (it == dims.nodeId2idx.end())
      return false;
    idx = it->second;
    return true;
  }
  inline bool valid_bounds(double lo, double hi)
  {
    return std::isfinite(lo) && std::isfinite(hi) && lo <= hi;
  }
  inline bool in_bounds(double x, double lo, double hi) { return x >= lo && x <= hi; }
} // namespace

bool VFMLoader::load_params(const XMLInput &dto, VFMState &state, std::string &err)
{
  state.params.clear();
  state.params.reserve(dto.parameters.size());

  for (const auto &p : dto.parameters)
  { // XMLInput::Param {name, init, lo, hi, scale}
    if (p.name.empty())
    {
      err = "Parameters: empty name";
      return false;
    }
    if (!std::isfinite(p.init))
    {
      err = "Parameters[" + p.name + "]: non-finite init";
      return false;
    }
    if (!valid_bounds(p.lo, p.hi))
    {
      err = "Parameters[" + p.name + "]: invalid bounds";
      return false;
    }
    if (!in_bounds(p.init, p.lo, p.hi))
    {
      err = "Parameters[" + p.name + "]: init out of bounds";
      return false;
    }
    if (!std::isfinite(p.scale) || p.scale == 0.0)
    {
      err = "Parameters[" + p.name + "]: invalid scale";
      return false;
    }

    VFMParam q;
    q.spec = VFMParamSpec{p.name, p.init, p.lo, p.hi, p.scale};
    q.value = p.init; // start at init
    state.params.push_back(q);
  }
  return true;
}

bool VFMLoader::load_measuredU(const XMLInput &dto, const MeshDims &dims, MeasuredData &out, std::string &err)
{
  out = MeasuredData{};
  out.setNodalSize(dims.nNodes);

  for (const auto &ts : dto.measuredU)
  {
    const TimeIdx t = out.addTime();
    for (const auto &s : ts.nodes)
    {
      std::size_t i;
      if (!map_node(dims, s.id, i))
      {
        err = "Unknown node id in measuredU: " + std::to_string(s.id);
        return false;
      }
      out.setU(t, i, vec3d(s.v.x, s.v.y, s.v.z));
    }
  }
  return true;
}

bool VFMLoader::load_virtualU(const XMLInput &dto, const MeshDims &dims, VirtualFields &out, std::string &err)
{
  out = VirtualFields{};
  out.resizeVF(dto.virtualU.size());
  out.setNodalSize(dims.nNodes);

  for (std::size_t v = 0; v < dto.virtualU.size(); ++v)
  {
    const auto &vf = dto.virtualU[v];
    for (const auto &ts : vf.times)
    {
      const TimeIdx t = out.addTime((VFIdx)v);
      for (const auto &s : ts.nodes)
      {
        std::size_t i;
        if (!map_node(dims, s.id, i))
        {
          err = "Unknown node id in virtualU: " + std::to_string(s.id);
          return false;
        }
        out.setU((VFIdx)v, t, i, vec3d(s.v.x, s.v.y, s.v.z));
      }
    }
  }
  return true;
}

bool VFMLoader::load_measuredF(const XMLInput &dto, const MeshDims &dims, MeasuredLoad &out, std::string &err)
{
  (void)dims;
  out = MeasuredLoad{};

  for (const auto &tl : dto.measuredLoads)
  {
    const TimeIdx t = out.addTime((double)tl.t);
    for (const auto &s : tl.loads)
    {
      out.addSurfaceLoad(t, s.surf, s.v);
    }
  }
  return true;
}
