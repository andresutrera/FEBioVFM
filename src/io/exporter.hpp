#pragma once
#include <string>

struct VFMProblem;

// Export computed displacement, deformation gradient, and stress fields to an XPLT file.
bool export_vfm_results(const VFMProblem& problem,
                        const std::string& filePath,
                        std::string& err);
