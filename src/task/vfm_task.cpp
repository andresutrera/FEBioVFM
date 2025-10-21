#include "vfm_task.h"
#include <FECore/log.h>
#include "optimization/vfm_solver.hpp"

VFMTask::VFMTask(FEModel *fem) : FECoreTask(fem) {}

bool VFMTask::Init(const char *xmlPath)
{
    m_inputPath = (xmlPath && *xmlPath) ? xmlPath : std::string{};
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

    if (!prepare_vfm_problem(*GetFEModel(), m_input, m_problem, err))
    {
        feLogError(err.c_str());
        return false;
    }

    feLog("Problem initialization complete.\n");
    return true;
}

bool VFMTask::Run()
{
    feLog("\n=== VFM RUN ===\n");
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
