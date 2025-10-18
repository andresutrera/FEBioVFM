// Minimal implementation scaffold for the Virtual Fields Method task. The goal
// is to keep enough structure so we can discuss how the dedicated VFMData.feb
// input will be consumed while leaving the heavy lifting for later steps.

#include "VFM.h"

#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <string>
#include "FEData.h"
#include "VFMKinematics.h"
#include "VFMExport.h"

/**
 * @brief Default constructor storing the FEBio model reference.
 *
 * @param fem FEBio model supplied by the kernel.
 * @note The context structure is left empty on purpose; it will be populated
 * once command line handling is expanded beyond simple file loading.
 */
FEVFMTask::FEVFMTask(FEModel* fem) : FECoreTask(fem), m_opt(fem)
{
}


/**
 * @brief Bootstrap the Virtual Fields Method task.
 *
 * The implementation announces the setup phase, loads optimization parameters
 * via FEOptimizeDataVFM, and initializes the registered variables. Any failure
 * encountered by the helper class is propagated upwards as a false return value.
 *
 * @param szfile Path to the optional VFM configuration file supplied by FEBio.
 * @return true if optimization data initialized successfully.
 * @note Additional context fields (e.g. measured data paths) are not yet parsed,
 * so @c m_context remains mostly empty after this call.
 */
bool FEVFMTask::Init(const char* szfile)
{
    feLog("V I R T U A L   F I E L D S   M E T H O D   (setup)\n");

	// read the data from the xml input file
	if (m_opt.Input(szfile) == false) return false;

	// do initialization
	bool ret = m_opt.Init();

	if (ret == false)
	{
		feLogErrorEx(m_opt.GetFEModel(), "Failed to initialize the optimization data.");
		return false;
	}

	std::string validationError;
	if (!VFMValidation::ValidateSolidDomains(*m_opt.GetFEModel(), validationError))
	{
		feLogErrorEx(m_opt.GetFEModel(), validationError.c_str());
		return false;
	}

	if (!VFMValidation::ValidateDisplacementCounts(*m_opt.GetFEModel(), m_opt, validationError))
	{
		feLogErrorEx(m_opt.GetFEModel(), validationError.c_str());
		return false;
	}

	std::string kinematicsError;
	if (!VFMKinematics::ComputeDeformationGradients(*m_opt.GetFEModel(),
		m_opt.MeasuredData(),
		m_opt.DeformationGradients(),
		kinematicsError))
	{
		feLogErrorEx(m_opt.GetFEModel(), kinematicsError.c_str());
		return false;
	}

	std::string plotPath;
	if (szfile && *szfile)
	{
		plotPath = szfile;
		std::size_t pos = plotPath.find_last_of('.');
		if (pos != std::string::npos)
			plotPath.replace(pos, std::string::npos, ".xplt");
		else
			plotPath.append(".xplt");
	}
	else
	{
		plotPath = "vfm_state.xplt";
	}

	std::string exportError;
	if (!ExportVFMKinematics(plotPath,
		*m_opt.GetFEModel(),
		m_opt.MeasuredData(),
		m_opt.VirtualData(),
		m_opt.DeformationGradients(),
		exportError))
	{
		feLogErrorEx(m_opt.GetFEModel(), exportError.c_str());
		return false;
	}
	else
	{
		feLog("VFM: exported kinematic snapshot to %s\n", plotPath.c_str());
	}

	return true;

}

/**
 * @brief Run the Virtual Fields Method task.
 *
 * @return true even though the body currently only emits TODO markers.
 * @note The success return keeps FEBio from aborting while the scaffolding is
 * developed, and the log output informs users about the missing pieces.
 */
bool FEVFMTask::Run()
{
    feLog("V I R T U A L   F I E L D S   M E T H O D   (run)\n");

    feLog("  TODO: parse measured displacements, virtual fields, options, and work history.\n");
    feLog("  TODO: assemble forward FEBio solve and apply Virtual Fields inversion.\n");

    return true;
}
