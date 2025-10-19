// Minimal implementation scaffold for the Virtual Fields Method task. The goal
// is to keep enough structure so we can discuss how the dedicated VFMData.feb
// input will be consumed while leaving the heavy lifting for later steps.

#include "VFM.h"

#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <string>
#include <vector>
#include "FEData.h"
#include "VFMKinematics.h"
#include "VFMExport.h"
#include "VFMStress.h"

namespace
{

void LogMatrix(FEModel& fem, const char* indent, const mat3d& m)
{
	feLogDebugEx(&fem, "%s[% .6e % .6e % .6e]", indent, m[0][0], m[0][1], m[0][2]);
	feLogDebugEx(&fem, "%s[% .6e % .6e % .6e]", indent, m[1][0], m[1][1], m[1][2]);
	feLogDebugEx(&fem, "%s[% .6e % .6e % .6e]", indent, m[2][0], m[2][1], m[2][2]);
}

void LogMatrix(FEModel& fem, const char* indent, const mat3ds& s)
{
	feLogDebugEx(&fem, "%s[% .6e % .6e % .6e]", indent, s.xx(), s.xy(), s.xz());
	feLogDebugEx(&fem, "%s[% .6e % .6e % .6e]", indent, s.xy(), s.yy(), s.yz());
	feLogDebugEx(&fem, "%s[% .6e % .6e % .6e]", indent, s.xz(), s.yz(), s.zz());
}

void LogParameterSummary(FEOptimizeDataVFM& opt)
{
	FEModel* fem = opt.GetFEModel();
	if (fem == nullptr) return;

	const int paramCount = opt.InputParameters();
	feLogDebugEx(fem, "  Parameters: %d", paramCount);
	if (paramCount == 0)
	{
		feLogDebugEx(fem, "    <none>");
		return;
	}

	for (int i = 0; i < paramCount; ++i)
	{
		FEInputParameterVFM* param = opt.GetInputParameter(i);
		if (param == nullptr) continue;

		feLogDebugEx(fem, "    %-20s init=%-12g min=%-12g max=%-12g",
			param->GetName().c_str(),
			param->InitValue(),
			param->MinValue(),
			param->MaxValue());
	}
}

void LogDisplacementHistory(FEModel& fem, const char* label, const DisplacementHistory& history)
{
	feLogDebugEx(&fem, "  %s displacements: %zu steps", label, history.Steps());
	if (history.Steps() == 0)
	{
		feLogDebugEx(&fem, "    <none>");
		return;
	}

	size_t stepIdx = 0;
	for (const auto& step : history)
	{
		const size_t nodeCount = step.displacements.Size();
		feLogDebugEx(&fem, "    [%02zu] t = %-12g nodes = %zu", stepIdx++, step.time, nodeCount);

		if (nodeCount == 0)
		{
			feLogDebugEx(&fem, "      <no displacement samples>");
			continue;
		}

		for (const NodeDisplacement& entry : step.displacements.Samples())
		{
			feLogDebugEx(&fem, "      node %6d : ux=%-12g uy=%-12g uz=%-12g",
				entry.id,
				entry.displacement[0],
				entry.displacement[1],
				entry.displacement[2]);
		}
	}
}

void LogLoadHistory(FEModel& fem, const MeasuredLoadHistory& history)
{
	feLogDebugEx(&fem, "  Measured loads: %zu steps", history.Steps());
	if (history.Steps() == 0)
	{
		feLogDebugEx(&fem, "    <none>");
		return;
	}

	size_t stepIdx = 0;
	for (const auto& step : history)
	{
		const auto& set = step.loads;
		feLogDebugEx(&fem, "    [%02zu] t = %-12g surfaces = %zu", stepIdx++, step.time, set.Size());

		if (set.Size() == 0)
		{
			feLogDebugEx(&fem, "      <no load samples>");
			continue;
		}

		for (const SurfaceLoadSample& sample : set.Samples())
		{
			feLogDebugEx(&fem, "      %-12s : Fx=%-12g Fy=%-12g Fz=%-12g",
				sample.id.c_str(),
				sample.load.x,
				sample.load.y,
				sample.load.z);
		}
	}
}

void LogVirtualFields(FEModel& fem, const VirtualDisplacementCollection& fields)
{
	feLogDebugEx(&fem, "  Virtual fields: %zu", fields.Size());
	if (fields.Empty())
	{
		feLogDebugEx(&fem, "    <none>");
		return;
	}

	for (size_t fieldIdx = 0; fieldIdx < fields.Size(); ++fieldIdx)
	{
		const auto& field = fields[fieldIdx];
		std::string label = "Virtual ";
		if (!field.id.empty())
		{
			label += "[" + field.id + "]";
		}
		else
		{
			label += "[#" + std::to_string(fieldIdx) + "]";
		}
		LogDisplacementHistory(fem, label.c_str(), field.history);
		++fieldIdx;
	}
}

void LogDeformationHistory(FEModel& fem, const DeformationGradientHistory& history)
{
	feLogDebugEx(&fem, "  Deformation gradients: %zu steps", history.Steps());
	if (history.Steps() == 0)
	{
		feLogDebugEx(&fem, "    <none>");
		return;
	}

	size_t stepIdx = 0;
	for (const auto& step : history)
	{
		const auto& elements = step.field.Data();
		feLogDebugEx(&fem, "    [%02zu] t = %-12g elements = %zu", stepIdx++, step.time, elements.size());

		if (elements.empty())
		{
			feLogDebugEx(&fem, "      <no deformation data>");
			continue;
		}

		for (const GaussPointDeformation& gp : elements)
		{
			feLogDebugEx(&fem, "      elem %6d : %zu gauss points", gp.elementId, gp.gradients.size());
			if (gp.gradients.empty())
			{
				feLogDebugEx(&fem, "        <no gradients>");
				continue;
			}

			size_t gpIdx = 0;
			for (const mat3d& F : gp.gradients)
			{
				feLogDebugEx(&fem, "        gp %02zu :", gpIdx);
				LogMatrix(fem, "          ", F);
				++gpIdx;
			}
		}
	}
}

void LogStressHistory(FEModel& fem, const StressHistory& history)
{
	feLogDebugEx(&fem, "  Stresses: %zu steps", history.Steps());
	if (history.Steps() == 0)
	{
		feLogDebugEx(&fem, "    <none>");
		return;
	}

	size_t stepIdx = 0;
	for (const auto& step : history)
	{
		const auto& elements = step.field.Data();
		feLogDebugEx(&fem, "    [%02zu] t = %-12g elements = %zu", stepIdx++, step.time, elements.size());

		if (elements.empty())
		{
			feLogDebugEx(&fem, "      <no stresses>");
			continue;
		}

		for (const GaussPointStress& gp : elements)
		{
			feLogDebugEx(&fem, "      elem %6d : %zu gauss points", gp.elementId, gp.stresses.size());
			if (gp.stresses.empty())
			{
				feLogDebugEx(&fem, "        <no stress tensors>");
				continue;
			}

			size_t gpIdx = 0;
			for (const mat3ds& sigma : gp.stresses)
			{
				feLogDebugEx(&fem, "        gp %02zu :", gpIdx);
				LogMatrix(fem, "          ", sigma);
				++gpIdx;
			}
		}
	}
}

void LogSetupDiagnostics(FEOptimizeDataVFM& opt)
{
	FEModel* fem = opt.GetFEModel();
	if (fem == nullptr) return;

	feLogDebugEx(fem, "---- VFM Diagnostics (setup) -------------------------");
	LogParameterSummary(opt);
	LogDisplacementHistory(*fem, "Measured", opt.MeasuredHistory());
	LogVirtualFields(*fem, opt.VirtualFields());
	LogLoadHistory(*fem, opt.MeasuredLoads());
	LogDeformationHistory(*fem, opt.DeformationHistory());
}

void LogStressDiagnostics(FEOptimizeDataVFM& opt)
{
	FEModel* fem = opt.GetFEModel();
	if (fem == nullptr) return;

	feLogDebugEx(fem, "---- VFM Diagnostics (stresses) ----------------------");
	LogStressHistory(*fem, opt.StressTimeline());
}

bool BuildStressHistory(FEOptimizeDataVFM& opt, std::string& errorMessage)
{
	auto& defHistory = opt.DeformationHistory();
	auto& stressHistory = opt.StressTimeline();

	stressHistory.Clear();
	stressHistory.Reserve(defHistory.Steps());

	if (defHistory.Empty()) return true;

	for (const auto& defStep : defHistory)
	{
		auto& stressStep = stressHistory.AddStep(defStep.time);

		if (!VFMStress::ComputeCauchyStress(*opt.GetFEModel(),
			defStep.field,
			stressStep.field,
			errorMessage))
		{
			return false;
		}

		// feLog("VFM: computed stresses for t = %g\n", defStep.time);
	}

	return true;
}

} // namespace

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


	// compute deformation gradients for each time step
	const auto& measuredHistory = m_opt.MeasuredHistory();
	auto& defHistory = m_opt.DeformationHistory();
	defHistory.Clear();
	defHistory.Reserve(measuredHistory.Steps());

	for (const auto& measStep : measuredHistory)
	{
		const double t = measStep.time;
		auto& defStep = defHistory.AddStep(t);
		defStep.field.Clear();

		std::string kinematicsError;
		if (!VFMKinematics::ComputeDeformationGradients(*m_opt.GetFEModel(),
			measStep.displacements,
			defStep.field,
			kinematicsError))
		{
			feLogErrorEx(m_opt.GetFEModel(), kinematicsError.c_str());
			return false;
		}

		feLog("VFM: computed deformation gradients for t = %g\n", t);
	}


	if (!VFMValidation::ValidateDisplacementCounts(*m_opt.GetFEModel(), m_opt, validationError))
	{
		feLogErrorEx(m_opt.GetFEModel(), validationError.c_str());
		return false;
	}

	if (!VFMValidation::ValidateMeasuredLoads(*m_opt.GetFEModel(), m_opt, validationError))
	{
		feLogErrorEx(m_opt.GetFEModel(), validationError.c_str());
		return false;
	}

	std::string stressError;
	if (!BuildStressHistory(m_opt, stressError))
	{
		feLogErrorEx(m_opt.GetFEModel(), stressError.c_str());
		return false;
	}

	LogSetupDiagnostics(m_opt);
	LogStressDiagnostics(m_opt);

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
		m_opt.MeasuredHistory(),
		m_opt.VirtualFields(),
		m_opt.DeformationHistory(),
		m_opt.StressTimeline(),
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

	std::string stressError;
	if (!BuildStressHistory(m_opt, stressError))
	{
		feLogErrorEx(m_opt.GetFEModel(), stressError.c_str());
		return false;
	}

	LogStressDiagnostics(m_opt);

	feLog("VFM: stress reconstruction completed for %zu time steps.\n", m_opt.StressTimeline().Steps());

	return true;
}
