#pragma once

#include <string>

#include <FECore/FECoreTask.h>
#include "FEData.h"

/**
 * @brief Lightweight context shared across the Virtual Fields Method workflow.
 *
 * The struct currently records the path to the secondary VFM data file and an
 * optional human readable string that can be surfaced in logs or UI. Additional
 * fields can be added incrementally as the plugin matures.
 */
struct FEVFMContext
{
    std::string dataFile;      ///< Path passed via "-task=\"VFM\" VFMData.feb".
    std::string description;   ///< Human readable summary of the current study.
};

/**
 * @brief FEBio task entry point that orchestrates the Virtual Fields Method.
 *
 * The task integrates with FEBio's command line via @c -task="VFM" <file>, reads
 * optimization parameters, and will ultimately run the forward solve followed by
 * the Virtual Fields inversion. The current implementation focuses on bootstrapping
 * the optimization container and leaves actual inversion steps as future work.
 *
 * @note Documentation intentionally calls out the placeholders so the generated
 * Doxygen pages help stakeholders understand the roadmap.
 */
class FEVFMTask : public FECoreTask
{
public:
    /**
     * @brief Create a task instance bound to the supplied FEBio model.
     * @param fem FEBio model provided by the hosting FEBio application.
     */
    explicit FEVFMTask(FEModel* fem);

    /**
     * @brief Initialize the task from the command line argument supplied to FEBio.
     * @param szfile Path to the optional VFM configuration file.
     * @return true if the optimization data was parsed successfully.
     *
     * @note At present the method delegates parsing to FEOptimizeDataVFM and only
     * logs minimal progress. Future revisions will populate @c m_context.
     */
    bool Init(const char* szfile) override;

    /**
     * @brief Execute the task after all initialization steps have completed.
     * @return true if the run completed without fatal errors.
     *
     * @note The current body simply prints TODO markers, making the unfinished
     * nature of the workflow explicit in the generated documentation.
     */
    bool Run() override;

    /**
     * @brief Fetch the cached Virtual Fields context structure.
     */
    const FEVFMContext& Context() const { return m_context; }

private:
    FEVFMContext m_context;      ///< Task-wide context populated during Init().
    FEOptimizeDataVFM	m_opt;    ///< Optimization data wrapper responsible for parsing input.
};

namespace FEBioVFM
{
}
