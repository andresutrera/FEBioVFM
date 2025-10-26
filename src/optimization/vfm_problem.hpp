#pragma once
#include <string>
#include "build/mesh_info.hpp"
#include "build/surface_info.hpp"
#include "state/vfm_state.hpp"
#include "io/xml_reader.hpp"

class FEModel;

struct VFMProblem
{
  FEModel *fem = nullptr;
  MeshDims dims;
  MeshConn conn;
  MeshQuad quad;
  SurfaceMap surfaces;
  VFMState state;
  XMLInput::Options solverOptions;

  void reset();
};

bool prepare_vfm_problem(FEModel &fem,
                         const XMLInput &input,
                         VFMProblem &problem,
                         std::string &err);
