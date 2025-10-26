#pragma once
#include <vector>
#include <string>
#include <array>

enum class VFMOptimizationMethod
{
    Levmar,
    ConstrainedLevmar
};

struct VFMOptimizationOptions
{
    static constexpr std::size_t optionCount = 5;

    VFMOptimizationMethod method = VFMOptimizationMethod::ConstrainedLevmar;
    std::array<double, optionCount> values{};
    std::array<bool, optionCount> overrides{};
    int maxIterations = -1;
};

class InternalWorkAssembler;
struct VFMProblem;

// Runs a simple Levenberg–Marquardt fit that minimizes
// ‖InternalWork(params) - externalWork‖_2.
// On success, 'params' is updated in-place with the optimized values.
// Returns false when the model evaluation fails; 'err' then carries the reason.
bool run_vfm_levmar(std::vector<double>& params,
                    InternalWorkAssembler& internal,
                    const std::vector<double>& externalWork,
                    const std::vector<double>& lowerBounds,
                    const std::vector<double>& upperBounds,
                    const VFMOptimizationOptions& options,
                    int itmax,
                    std::string& err);

bool solve_vfm_problem(VFMProblem& problem,
                       int itmax,
                       std::string& err);
