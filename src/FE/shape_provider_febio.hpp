#pragma once
#include <vector>
#include "services/shape_provider.hpp"
#include "build/mesh_info.hpp"  // MeshConn, ElemRef
class FEBioShapeProvider final : public IShapeProvider {
public:
  explicit FEBioShapeProvider(const MeshConn& c): _c(c) {}
  const std::vector<size_t>& elemNodes(size_t e) const override { return _c.elemNodes[e]; }
  void gradN(size_t e, size_t g, std::vector<vec3d>& dNdx0) const override;
private:
  const MeshConn& _c;
};
