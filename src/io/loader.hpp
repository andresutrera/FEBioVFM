#pragma once
#include <string>
#include "io/xml_reader.hpp"          // XMLInput (int times, ids)
#include "build/mesh_info.hpp"        // MeshDims
#include "domain/vfm_displacements.hpp" // MeasuredData, VirtualFields, MeasuredLoad


struct XMLInput;
struct MeshDims;
struct MeasuredData;
struct VirtualFields;
struct MeasuredLoad;
struct VFMState;


namespace VFMLoader {
  bool load_params(const XMLInput& dto, VFMState& state, std::string& err);
  bool load_measuredU (const XMLInput& dto, const MeshDims& dims, MeasuredData& out, std::string& err);
  bool load_virtualU  (const XMLInput& dto, const MeshDims& dims, VirtualFields& out, std::string& err);
  bool load_measuredF (const XMLInput& dto, const MeshDims& dims, MeasuredLoad& out, std::string& err);
}
