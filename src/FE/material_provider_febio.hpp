// fe/material_provider_febio.hpp
#pragma once
#include <string>
#include "services/stress_eval.hpp"
#include "build/mesh_info.hpp"
class FEBioMaterialProvider final : public IMaterialProvider {
public:
  FEBioMaterialProvider(const MeshConn& c) : _c(c) {}
  bool evalCauchy(size_t e, size_t g, const mat3d& F, mat3ds& sigma, std::string& err) const override;
private:
  const MeshConn& _c; // has ElemRef {FESolidDomain*, local}
};
