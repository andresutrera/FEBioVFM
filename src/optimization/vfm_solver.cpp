#include "optimization/vfm_solver.hpp"
#include "vfm/InternalWork.hpp"
#include "levmar/levmar.h"

namespace
{
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
