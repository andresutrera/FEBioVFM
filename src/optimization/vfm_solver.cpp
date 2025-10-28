#include "optimization/vfm_solver.hpp"
#include "optimization/vfm_problem.hpp"
#include "vfm/InternalWork.hpp"
#include "vfm/ExternalVirtualWork.hpp"
#include "FE/params_febio.hpp"
#include "FE/material_provider_febio.hpp"
#include "services/stress_eval.hpp"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <utility>
#include "levmar/levmar.h"
#include "diag/felog_bridge.hpp"
#include "diag/printers/param_table.hpp"
#include <array>
#include <atomic>
#include <csignal>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace
{
    constexpr int kDefaultMaxIterations = 100;

    std::atomic_bool g_levmarInterrupt{false};

    void handle_sigint(int)
    {
        g_levmarInterrupt.store(true);
    }

    struct ScopedSigintHandler
    {
        using Handler = void (*)(int);
        Handler previous = nullptr;
        bool installed = false;

        ScopedSigintHandler()
        {
#ifdef SIGINT
            previous = std::signal(SIGINT, handle_sigint);
            installed = true;
#endif
            g_levmarInterrupt.store(false);
        }

        ~ScopedSigintHandler()
        {
#ifdef SIGINT
            if (installed)
                std::signal(SIGINT, previous);
#endif
        }
    };

    struct LevmarContext
    {
        InternalWorkAssembler *internal = nullptr;
        std::string *err = nullptr;
        bool failed = false;
        const std::vector<double> *external = nullptr;
        bool logEvaluations = false;
        std::size_t evalCount = 0;
        std::atomic_bool *cancelFlag = nullptr;
    };

    void lm_internal_eval(double *p, double *hx, int m, int n, void *adata)
    {
        auto *ctx = static_cast<LevmarContext *>(adata);
        if (!ctx || !ctx->internal)
        {
            for (int i = 0; i < n; ++i)
                hx[i] = 0.0;
            if (ctx && ctx->err)
                *ctx->err = "invalid levmar context";
            if (ctx)
                ctx->failed = true;
            return;
        }

        if (ctx->failed)
        {
            for (int i = 0; i < n; ++i)
                hx[i] = 0.0;
            return;
        }

        if (ctx->cancelFlag && ctx->cancelFlag->load())
        {
            if (ctx->err && ctx->err->empty())
                *ctx->err = "optimization interrupted";
            ctx->failed = true;
            for (int i = 0; i < n; ++i)
                hx[i] = 0.0;
            return;
        }

        std::vector<double> params(p, p + m);
        std::string localErr;
        std::vector<double> iw = (*ctx->internal)(params, localErr);
        if (iw.size() != static_cast<std::size_t>(n))
        {
            if (ctx->err)
                *ctx->err = localErr.empty() ? "internal work dimension mismatch" : localErr;
            ctx->failed = true;
            for (int i = 0; i < n; ++i)
                hx[i] = 0.0;
            return;
        }

        if (!localErr.empty() && ctx->err)
            *ctx->err = localErr;

        if (ctx->logEvaluations)
        {
            double cost = 0.0;
            bool costValid = ctx->external && ctx->external->size() == static_cast<std::size_t>(n);
            if (costValid)
            {
                for (int i = 0; i < n; ++i)
                {
                    const double diff = iw[static_cast<std::size_t>(i)] - ctx->external->at(static_cast<std::size_t>(i));
                    cost += diff * diff;
                }
                cost *= 0.5;
            }

            std::ostringstream oss;
            oss.setf(std::ios::scientific);
            oss << std::setprecision(6);
            oss << "LM eval " << (++ctx->evalCount);
            if (costValid)
                oss << " | cost=" << cost;
            else
                oss << " | cost=N/A";
            oss << " | params=[";
            for (int i = 0; i < m; ++i)
            {
                if (i > 0)
                    oss << ", ";
                oss << params[static_cast<std::size_t>(i)];
            }
            oss << "]";
            feLog("%s\n", oss.str().c_str());
        }

        for (int i = 0; i < n; ++i)
            hx[i] = iw[static_cast<std::size_t>(i)];
    }

    VFMOptimizationOptions make_solver_options(const XMLInput::Options &src)
    {
        VFMOptimizationOptions opts;

        if (src.present && src.type == XMLInput::Options::Type::Levmar)
            opts.method = VFMOptimizationMethod::Levmar;

        if (src.tau.set)
        {
            opts.overrides[0] = true;
            opts.values[0] = src.tau.value;
        }
        if (src.gradTol.set)
        {
            opts.overrides[1] = true;
            opts.values[1] = src.gradTol.value;
        }
        if (src.stepTol.set)
        {
            opts.overrides[2] = true;
            opts.values[2] = src.stepTol.value;
        }
        if (src.objTol.set)
        {
            opts.overrides[3] = true;
            opts.values[3] = src.objTol.value;
        }
        if (src.diffScale.set)
        {
            opts.overrides[4] = true;
            opts.values[4] = src.diffScale.value;
        }
        opts.maxIterations = kDefaultMaxIterations;
        if (src.maxIters.set && src.maxIters.value > 0.0)
            opts.maxIterations = static_cast<int>(src.maxIters.value);

        return opts;
    }
} // namespace

bool run_vfm_levmar(std::vector<double> &params,
                    InternalWorkAssembler &internal,
                    const std::vector<double> &externalWork,
                    const std::vector<double> &lowerBounds,
                    const std::vector<double> &upperBounds,
                    const VFMOptimizationOptions &options,
                    std::string &err)
{
    err.clear();

    if (params.empty())
        return true;
    if (externalWork.empty())
    {
        err = "external work vector is empty";
        return false;
    }

    const int m = static_cast<int>(params.size());
    const int n = static_cast<int>(externalWork.size());

    const bool useBounds = options.method == VFMOptimizationMethod::ConstrainedLevmar;
    if (useBounds)
    {
        if (lowerBounds.size() != params.size() || upperBounds.size() != params.size())
        {
            err = "bounds size mismatch";
            return false;
        }
    }

    std::vector<double> x = externalWork;

    LevmarContext ctx;
    ctx.internal = &internal;
    ctx.err = &err;
    ctx.external = &x;
    ctx.logEvaluations = true;
    ctx.evalCount = 0;
    ctx.cancelFlag = &g_levmarInterrupt;

    ScopedSigintHandler sigintGuard;

    std::array<double, VFMOptimizationOptions::optionCount> opts = {1e-3, 1e-12, 1e-12, 1e-15, -1.0};
    for (std::size_t i = 0; i < opts.size(); ++i)
    {
        if (options.overrides[i])
            opts[i] = options.values[i];
    }

    std::vector<double> lb;
    std::vector<double> ub;
    if (useBounds)
    {
        lb.assign(lowerBounds.begin(), lowerBounds.end());
        ub.assign(upperBounds.begin(), upperBounds.end());
    }

    const std::size_t workSize = useBounds ? static_cast<std::size_t>(LM_BC_DIF_WORKSZ(m, n))
                                           : static_cast<std::size_t>(LM_DIF_WORKSZ(m, n));
    std::vector<double> work(workSize > 0 ? workSize : 1u, 0.0);
    double info[LM_INFO_SZ] = {};
    ctx.failed = false;

    int status;
    const int maxIterations = options.maxIterations > 0 ? options.maxIterations : kDefaultMaxIterations;
    if (useBounds)
    {
        status = dlevmar_bc_dif(&lm_internal_eval,
                                params.data(), x.data(),
                                m, n,
                                lb.data(), ub.data(), nullptr,
                                maxIterations,
                                opts.data(), info,
                                work.data(), nullptr,
                                &ctx);
    }
    else
    {
        status = dlevmar_dif(&lm_internal_eval,
                             params.data(), x.data(),
                             m, n,
                             maxIterations,
                             opts.data(), info,
                             work.data(), nullptr,
                             &ctx);
    }

    if (ctx.failed)
    {
        if (err.empty() && g_levmarInterrupt.load())
            err = "optimization interrupted";
        return false;
    }
    if (status < 0)
    {
        if (err.empty())
            err = "levmar failed";
        return false;
    }

    feLog("\nLEV-MAR SUMMARY\n");
    feLog("  Initial cost  : %.6e\n", info[0]);
    feLog("  Final cost    : %.6e\n", info[1]);
    feLog("  ||J^T e||_inf : %.6e\n", info[2]);
    feLog("  ||dx||        : %.6e\n", info[3]);
    feLog("  mu/max diag   : %.6e\n", info[4]);
    feLog("  Iterations    : %d\n", static_cast<int>(info[5]));
    feLog("  Stop reason   : %d\n", static_cast<int>(info[6]));
    feLog("  Function evals: %d\n", static_cast<int>(info[7]));
    feLog("  Jacobians     : %d\n", static_cast<int>(info[8]));
    feLog("  Linear solves : %d\n", static_cast<int>(info[9]));

    return true;
}

bool solve_vfm_problem(VFMProblem &problem,
                       std::string &err)
{
    diag::ScopedFEBind bind(problem.fem);
    err.clear();
    if (problem.fem == nullptr)
    {
        err = "VFM problem not initialized";
        return false;
    }

    if (problem.state.params.empty())
    {
        feLog("No parameters to optimize.\n");
        return true;
    }

    FEModel &fem = *problem.fem;
    FEBioParameterApplier paramApplier(fem, problem.state);
    FEBioMaterialProvider mat(problem.conn);

    auto paramSetter = [&](const std::vector<double> &values, std::string &perr) -> bool
    {
        return paramApplier.apply(values, perr);
    };

    auto computeStress = [&](std::string &perr) -> bool
    {
        if (!StressEval::cauchy(problem.state.def, problem.state.stresses, mat, perr))
            return false;
        return StressEval::first_piola(problem.state.def, problem.state.stresses, problem.state.stresses, perr);
    };

    auto toVirtGrad = [](const mat3d &Fstar)
    {
        mat3d G = Fstar;
        G[0][0] -= 1.0;
        G[1][1] -= 1.0;
        G[2][2] -= 1.0;
        return G;
    };

    InternalWorkAssembler internal(problem.dims, problem.quad,
                                   problem.state.vdef, problem.state.stresses,
                                   paramSetter, computeStress, toVirtGrad);

    ExternalVirtualWorkAssembler external(problem.surfaces,
                                          problem.state.virtuals,
                                          problem.state.loads);

    std::vector<double> ew = external(err);
    if (!err.empty())
        return false;
    if (ew.empty())
    {
        feLog("External work vector empty. Nothing to optimize.\n");
        return true;
    }

    std::vector<double> params(problem.state.params.size());
    for (std::size_t i = 0; i < params.size(); ++i)
        params[i] = problem.state.params[i].value;

    std::vector<double> lower(params.size());
    std::vector<double> upper(params.size());
    for (std::size_t i = 0; i < params.size(); ++i)
    {
        lower[i] = problem.state.params[i].spec.lo;
        upper[i] = problem.state.params[i].spec.hi;
    }

    VFMOptimizationOptions solverOpts = make_solver_options(problem.solverOptions);

    if (!run_vfm_levmar(params, internal, ew, lower, upper, solverOpts, err))
        return false;

    if (!paramApplier.apply(params, err))
        return false;

    for (std::size_t i = 0; i < problem.state.params.size(); ++i)
        problem.state.params[i].value = params[i];

    if (!computeStress(err))
        return false;

    auto noopSetter = [](const std::vector<double> &, std::string &perr) -> bool
    {
        perr.clear();
        return true;
    };
    auto noopStress = [](std::string &perr) -> bool
    {
        perr.clear();
        return true;
    };

    InternalWorkAssembler finalInternal(problem.dims, problem.quad,
                                        problem.state.vdef, problem.state.stresses,
                                        noopSetter, noopStress, toVirtGrad);

    std::string iwErr;
    std::vector<double> iw = finalInternal(params, iwErr);
    if (!iwErr.empty())
    {
        err = iwErr;
        return false;
    }

    feLog("\nWork arrays after optimization:\n");
    feLog("External virtual work (size=%zu):\n", ew.size());
    for (std::size_t i = 0; i < ew.size(); ++i)
        feLog("  evw[%zu] = %.6e\n", i, ew[i]);
    feLog("Internal virtual work (size=%zu):\n", iw.size());
    for (std::size_t i = 0; i < iw.size(); ++i)
        feLog("  ivw[%zu] = %.6e\n", i, iw[i]);

    if (problem.solverOptions.saveVirtualWorkSet && !problem.solverOptions.saveVirtualWork.empty())
    {
        const std::size_t vfCount = problem.state.virtuals.nVF();
        const std::size_t iwSize = iw.size();
        const std::size_t ewSize = ew.size();

        if (vfCount == 0 || iwSize == 0 || ewSize == 0)
        {
            feLogWarning("Virtual work export requested, but no data available.\n");
        }
        else if (iwSize != ewSize)
        {
            err = "internal/external virtual work size mismatch.";
            return false;
        }
        else if (iwSize % vfCount != 0)
        {
            err = "virtual work data size is not divisible by the number of virtual fields.";
            return false;
        }
        else
        {
            const std::size_t timeCount = iwSize / vfCount;
            const std::size_t loadTimes = problem.state.loads.nTimes();
            const bool haveLoadTimes = loadTimes == timeCount && timeCount > 0;

            std::ofstream out(problem.solverOptions.saveVirtualWork);
            if (!out)
            {
                err = "failed to open virtual work output file: " + problem.solverOptions.saveVirtualWork;
                return false;
            }

            out.setf(std::ios::scientific);
            out << std::setprecision(6);

            out << "#Step";
            for (std::size_t v = 0; v < vfCount; ++v)
                out << ", IVW" << (v + 1);
            for (std::size_t v = 0; v < vfCount; ++v)
                out << ", EVW" << (v + 1);
            out << '\n';

            for (std::size_t t = 0; t < timeCount; ++t)
            {
                if (haveLoadTimes)
                {
                    const auto &frame = problem.state.loads.frame(static_cast<TimeIdx>(t));
                    out << frame.time;
                }
                else
                {
                    out << "t" << t;
                }

                for (std::size_t v = 0; v < vfCount; ++v)
                {
                    const std::size_t idx = v * timeCount + t;
                    out << ", " << iw[idx];
                }
                for (std::size_t v = 0; v < vfCount; ++v)
                {
                    const std::size_t idx = v * timeCount + t;
                    out << ", " << ew[idx];
                }
                out << '\n';
            }

            if (!out.good())
            {
                err = "failed while writing virtual work output file: " + problem.solverOptions.saveVirtualWork;
                return false;
            }

            feLog("Virtual work saved to %s\n", problem.solverOptions.saveVirtualWork.c_str());
        }
    }

    feLog("\n");
    diag::printers::ParameterTable(problem.state.params, "FINAL PARAMETERS", 6);
    feLog("\n");
    return true;
}
