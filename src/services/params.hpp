#pragma once
#include <vector>
#include <string>
struct IParameterApplier {
  virtual ~IParameterApplier() = default;
  // p are unscaled decision variables in optimizer space
  virtual bool apply(const std::vector<double>& p, std::string& err) = 0;
};