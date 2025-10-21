#pragma once
#include <string>
#include <vector>
#include "build/mesh_info.hpp"            // MeshConn, MeshQuad
#include "domain/vfm_displacements.hpp"   // MeasuredData, VirtualFields
#include "domain/vfm_tensors.hpp"   // Deformations, VirtualDeformations
#include "services/shape_provider.hpp"


namespace Kinematics {
    // Assemble F(t,e,g) from measured nodal u
    bool compute_measured (const MeshQuad& quad, const IShapeProvider& shp,
                           const MeasuredData& U, Deformations& Fout,
                           bool checkDet, std::string& err);

    // Assemble F_v(v,t,e,g) from virtual nodal u
    bool compute_virtuals (const MeshQuad& quad, const IShapeProvider& shp,
                           const VirtualFields& UV, VirtualDeformations& FVout,
                           bool checkDet, std::string& err);
}
