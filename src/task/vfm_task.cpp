#include "vfm_task.h"
#include <FECore/log.h>
#include "optimization/vfm_solver.hpp"
#include "diag/felog_bridge.hpp"
#include "diag/printers/param_table.hpp"

VFMTask::VFMTask(FEModel *fem) : FECoreTask(fem) {}

bool VFMTask::Init(const char *xmlPath)
{
    m_inputPath = (xmlPath && *xmlPath) ? xmlPath : std::string{};
    diag::ScopedFEBind bind(GetFEModel());
    feLog("\n");
    feLog("===========================================================================\n");
    feLog("                        VIRTUAL FIELDS METHOD (VFM)                        \n");
    feLog("===========================================================================\n");
    feLog("\n");
    feLog("...........................................................................\n");
    feLog("                                   SETUP                                   \n");
    feLog("...........................................................................\n");
    feLog("\n\n");

    std::string err;
    VFMXmlReader rdr;
    if (!rdr.read(m_inputPath.c_str(), m_input, err))
    {
        feLogError(err.c_str());
        return false;
    }
    feLog("Success reading input files.\n");

    if (!prepare_vfm_problem(*GetFEModel(), m_input, m_problem, err))
    {
        feLogError(err.c_str());
        return false;
    }

    feLog("Problem initialization complete.\n");

    feLog("");
    diag::printers::ParameterTable(m_problem.state.params, "INITIAL PARAMETERS", 6);

    return true;
}

bool VFMTask::Run()
{
    feLog("...........................................................................\n");
    feLog("                                    RUN                                    \n");
    feLog("...........................................................................\n");
    feLog("\n");
    feLog("\n");

    std::string err;
    if (m_problem.fem == nullptr)
    {
        feLogError("VFM problem not initialized.");
        return false;
    }
    if (!solve_vfm_problem(m_problem, 100, err))
    {
        if (!err.empty())
            feLogError(err.c_str());
        else
            feLogError("VFM solver failed.");
        return false;
    }

    return true;
}
