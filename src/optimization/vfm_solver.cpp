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

namespace
{
struct ProblemLogger
{
    explicit ProblemLogger(FEModel &model) : fem(&model) {}
    FEModel *GetFEModel() const { return fem; }

    template <typename... Args>
    void info(const char *fmt, Args &&...args) const
    {
        feLog(fmt, std::forward<Args>(args)...);
    }

private:
    FEModel *fem = nullptr;
};

struct LevmarContext
{
    InternalWorkAssembler *internal = nullptr;
    std::string *err = nullptr;
    bool failed = false;
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

    for (int i = 0; i < n; ++i)
        hx[i] = iw[static_cast<std::size_t>(i)];
}
} // namespace

bool run_vfm_levmar(std::vector<double> &params,
                    InternalWorkAssembler &internal,
                    const std::vector<double> &externalWork,
                    int itmax,
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

    std::vector<double> x = externalWork;

    LevmarContext ctx;
    ctx.internal = &internal;
    ctx.err = &err;

    double opts[LM_OPTS_SZ];
    opts[0] = 1e-3;  // tau
    opts[1] = 1e-12; // eps1
    opts[2] = 1e-12; // eps2
    opts[3] = 1e-15; // eps3
    opts[4] = -1.0;  // delta (<0 → use internal)

    double info[LM_INFO_SZ];
    dlevmar_dif(&lm_internal_eval,
                params.data(), x.data(),
                m, n,
                itmax,
                opts, info,
                nullptr, nullptr,
                &ctx);

    if (ctx.failed)
        return false;

    return true;
}

bool solve_vfm_problem(VFMProblem &problem,
                       int itmax,
                       std::string &err)
{
    err.clear();
    if (problem.fem == nullptr)
    {
        err = "VFM problem not initialized";
        return false;
    }

    ProblemLogger logger(*problem.fem);

    if (problem.state.params.empty())
    {
        logger.info("No parameters to optimize.\n");
        return true;
    }

    FEModel &fem = *problem.fem;
    FEBioParameterApplier paramApplier(fem, problem.state);
    FEBioMaterialProvider mat(problem.conn);

    auto paramSetter = [&](const std::vector<double> &values, std::string &perr) -> bool {
        return paramApplier.apply(values, perr);
    };

    auto computeStress = [&](std::string &perr) -> bool {
        if (!StressEval::cauchy(problem.state.def, problem.state.stresses, mat, perr))
            return false;
        return StressEval::first_piola(problem.state.def, problem.state.stresses, problem.state.stresses, perr);
    };

    auto toVirtGrad = [](const mat3d &Fstar) {
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
        logger.info("External work vector empty. Nothing to optimize.\n");
        return true;
    }

    std::vector<double> params(problem.state.params.size());
    for (std::size_t i = 0; i < params.size(); ++i)
        params[i] = problem.state.params[i].value;

    if (!run_vfm_levmar(params, internal, ew, itmax, err))
        return false;

    if (!paramApplier.apply(params, err))
        return false;

    for (std::size_t i = 0; i < problem.state.params.size(); ++i)
        problem.state.params[i].value = params[i];

    if (!computeStress(err))
        return false;

    logger.info("Optimized parameters:\n");
    for (const auto &q : problem.state.params)
        logger.info("  %s = %.6g\n", q.spec.name.c_str(), q.value);

    logger.info("Optimization complete.\n");
    return true;
}
