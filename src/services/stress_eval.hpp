// services/stress_eval.hpp
#pragma once
#include <string>
#include "domain/vfm_tensors.hpp"

struct IMaterialProvider {
  virtual ~IMaterialProvider() = default;
  // Prepare a scratch material point for element e, gp g, set F, return Cauchy σ
  virtual bool evalCauchy(size_t e, size_t g, const mat3d& F, mat3ds& sigma, std::string& err) const = 0;
};

namespace StressEval {
  // σ(t,e,g) from F(t,e,g)
  bool cauchy(const Deformations& F, Stresses& out, const IMaterialProvider& mat, std::string& err);

  // P = J σ F^{-T}
  bool first_piola(const Deformations& F, const Stresses& S, Stresses& outP, std::string& err);
}
