// src/services/shape_provider.hpp
#pragma once
#include <vector>
#include <cstddef>
#include <FECore/vec3d.h> // or your vec3d
struct IShapeProvider {
  virtual ~IShapeProvider() = default;
  virtual const std::vector<size_t>& elemNodes(size_t e) const = 0;
  virtual void gradN(size_t e, size_t g, std::vector<vec3d>& dNdx0) const = 0;
};
