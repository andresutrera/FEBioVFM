// diag/felog_bridge.hpp
#pragma once
#include <FECore/FEModel.h>
#include <FECore/log.h>

namespace diag
{
    inline thread_local FEModel *tls_fem = nullptr;
    inline void set_current_fem(FEModel *m) { tls_fem = m; }

    // optional RAII
    struct ScopedFEBind
    {
        FEModel *prev;
        explicit ScopedFEBind(FEModel *m) : prev(tls_fem) { tls_fem = m; }
        ~ScopedFEBind() { tls_fem = prev; }
    };
} // namespace diag

// Expose exactly what feLog’s macro expects
inline FEModel *GetFEModel() { return diag::tls_fem; }
