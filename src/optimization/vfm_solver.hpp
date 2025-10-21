#pragma once
#include <vector>
#include <string>

class InternalWorkAssembler;
struct VFMProblem;

// Runs a simple Levenberg–Marquardt fit that minimizes
// ‖InternalWork(params) - externalWork‖_2.
// On success, 'params' is updated in-place with the optimized values.
// Returns false when the model evaluation fails; 'err' then carries the reason.
bool run_vfm_levmar(std::vector<double>& params,
                    InternalWorkAssembler& internal,
                    const std::vector<double>& externalWork,
                    int itmax,
                    std::string& err);

bool solve_vfm_problem(VFMProblem& problem,
                       int itmax,
                       std::string& err);
