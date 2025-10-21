#pragma once

#include <string>

#include <FECore/FECoreTask.h>
#include "FEData.h"
#include "VFMValidation.h"

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
class VFMTask : public FECoreTask
{
public:
    /**
     * @brief Create a task instance bound to the supplied FEBio model.
     * @param fem FEBio model provided by the hosting FEBio application.
     */
    explicit VFMTask(FEModel *fem);

    /**
     * @brief Initialize the task from the command line argument supplied to FEBio.
     * @param szfile Path to the optional VFM configuration file.
     * @return true if the optimization data was parsed successfully.
     *
     * @note At present the method delegates parsing to FEOptimizeDataVFM and only
     * logs minimal progress. Future revisions will populate @c m_context.
     */
    bool Init(const char *szfile) override;

    /**
     * @brief Execute the task after all initialization steps have completed.
     * @return true if the run completed without fatal errors.
     *
     * @note The current body simply prints TODO markers, making the unfinished
     * nature of the workflow explicit in the generated documentation.
     */
    bool Run() override;

private:
    bool LoadInput(const char *szfile);
    bool InitializeParameters();
    bool ValidateFEModel();
    bool ComputeMeasuredKinematics();
    bool ComputeVirtualKinematics();
    bool ValidateDataConsistency();
    bool ComputeExternalVirtualWork();
    bool LogDiagnostics();
    bool ExportState(const char *szfile);

    void LogVector(const char *tag, const std::vector<double> &v);
    void LogParameterTable(std::vector<FEInputParameterVFM *> &vars, const char *title, int precision);

    FEOptimizeDataVFM m_opt; ///< Optimization data wrapper responsible for parsing input.
    std::string m_inputFile; ///< Original VFM data path used to derive export filenames.
};

namespace FEBioVFM
{
}
