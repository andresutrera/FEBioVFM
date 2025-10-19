#include <unordered_map>
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
#include <FECore/FEFacetSet.h>

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

void LogVirtualDeformationHistories(FEModel& fem, const VirtualDeformationGradientCollection& fields)
{
	feLogDebugEx(&fem, "  Virtual deformation gradients: %zu", fields.Size());
	if (fields.Empty())
	{
		feLogDebugEx(&fem, "    <none>");
		return;
	}

	size_t idx = 0;
	for (const auto& field : fields.Data())
	{
		std::string label = "Virtual Def ";
		if (!field.id.empty())
		{
			label += "[" + field.id + "]";
		}
		else
		{
			label += "[#" + std::to_string(idx) + "]";
		}
		feLogDebugEx(&fem, label.c_str());
		LogDeformationHistory(fem, field.history);
		++idx;
	}
}

void LogVirtualExternalWork(FEModel& fem, const FEOptimizeDataVFM& opt)
{
	const auto& workHistories = opt.VirtualExternalWork();
	feLogDebugEx(&fem, "  Virtual external work: %zu fields", workHistories.size());
	if (workHistories.empty())
	{
		feLogDebugEx(&fem, "    <none>");
		return;
	}

	std::vector<double> times;
	for (const auto& step : opt.MeasuredLoads().StepsRef())
	{
		times.push_back(step.time);
	}

	size_t fieldIdx = 0;
	const auto& virtualFields = opt.VirtualFields();
	for (const auto& history : workHistories)
	{
		std::string label = "Virtual Work ";
		if (fieldIdx < virtualFields.Size() && !virtualFields[fieldIdx].id.empty())
		{
			label += "[" + virtualFields[fieldIdx].id + "]";
		}
		else
		{
			label += "[#" + std::to_string(fieldIdx) + "]";
		}

		feLogDebugEx(&fem, "    %s", label.c_str());
		const auto& work = history.work;
		for (size_t ti = 0; ti < work.size(); ++ti)
		{
			double t = (ti < times.size() ? times[ti] : static_cast<double>(ti));
			feLogDebugEx(&fem, "      t = %-12g Wext = %-12g", t, work[ti]);
		}
		++fieldIdx;
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
	LogVirtualDeformationHistories(*fem, opt.VirtualDeformationGradients());
	LogVirtualExternalWork(*fem, opt);
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

	if (!LoadInput(szfile)) return false;
	if (!InitializeOptimization()) return false;
	if (!ValidateModel()) return false;
	if (!ComputeMeasuredKinematics()) return false;
	if (!ComputeVirtualKinematics()) return false;
	if (!ValidateDataConsistency()) return false;
	if (!BuildStressHistoryStage()) return false;
	if (!ComputeExternalWork()) return false;
	LogDiagnostics();
	if (!ExportState(szfile)) return false;
	return true;
}

bool FEVFMTask::LoadInput(const char* szfile)
{
	if (m_opt.Input(szfile) == false) return false;
	return true;
}

bool FEVFMTask::InitializeOptimization()
{
	if (m_opt.Init() == false)
	{
		feLogErrorEx(m_opt.GetFEModel(), "Failed to initialize the optimization data.");
		return false;
	}
	return true;
}

bool FEVFMTask::ValidateModel()
{
	std::string validationError;
	if (!VFMValidation::ValidateSolidDomains(*m_opt.GetFEModel(), validationError))
	{
		feLogErrorEx(m_opt.GetFEModel(), validationError.c_str());
		return false;
	}
	return true;
}

bool FEVFMTask::ComputeMeasuredKinematics()
{
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
	return true;
}

bool FEVFMTask::ComputeVirtualKinematics()
{
	auto& virtualGradients = m_opt.VirtualDeformationGradients();
	virtualGradients.Clear();
	const auto& virtualFields = m_opt.VirtualFields();
	size_t fieldIdx = 0;
	for (const auto& field : virtualFields)
	{
		auto& outField = virtualGradients.Add(field.id);
		outField.history.Clear();
		outField.history.Reserve(field.history.Steps());

		for (const auto& step : field.history.StepsRef())
		{
			auto& gradStep = outField.history.AddStep(step.time);
			gradStep.field.Clear();

			std::string virtualError;
			if (!VFMKinematics::ComputeDeformationGradients(*m_opt.GetFEModel(),
				step.displacements,
				gradStep.field,
				virtualError))
			{
				const std::string name = field.id.empty() ? ("#" + std::to_string(fieldIdx)) : field.id;
				feLogErrorEx(m_opt.GetFEModel(), "Failed to compute virtual deformation gradients for field '%s' at t = %g: %s", name.c_str(), step.time, virtualError.c_str());
				return false;
			}
		}

		++fieldIdx;
	}
	return true;
}

bool FEVFMTask::ValidateDataConsistency()
{
	std::string validationError;
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

	return true;
}

bool FEVFMTask::BuildStressHistoryStage()
{
	std::string stressError;
	if (!BuildStressHistory(m_opt, stressError))
	{
		feLogErrorEx(m_opt.GetFEModel(), stressError.c_str());
		return false;
	}
	return true;
}

bool FEVFMTask::ComputeExternalWork()
{
	const auto& loadsHistory = m_opt.MeasuredLoads();
	const auto& virtualFields = m_opt.VirtualFields();
	auto& workHistories = m_opt.VirtualExternalWork();
	workHistories.clear();

	const auto& loadSteps = loadsHistory.StepsRef();
	const size_t T = loadSteps.size();
	const size_t N = virtualFields.Size();
	if ((T == 0) || (N == 0)) return true;

	std::vector<double> times(T, 0.0);
	for (size_t t = 0; t < T; ++t)
	{
		times[t] = loadSteps[t].time;
	}

	std::vector<std::string> surfaces;
	for (size_t t = 0; t < T; ++t)
	{
		const auto& samples = loadSteps[t].loads.Samples();
		for (size_t s = 0; s < samples.size(); ++s)
		{
			const std::string& name = samples[s].id;
			bool exists = false;
			for (size_t k = 0; k < surfaces.size(); ++k)
			{
				if (surfaces[k] == name)
				{
					exists = true;
					break;
				}
			}
			if (!exists) surfaces.push_back(name);
		}
	}

	const size_t K = surfaces.size();
	std::vector<std::vector<int>> surfaceNodes(K);
	FEMesh& mesh = m_opt.GetFEModel()->GetMesh();
	for (size_t k = 0; k < K; ++k)
	{
		const std::string& sid = surfaces[k];
		std::vector<int>& nodes = surfaceNodes[k];

		FESurface* surface = mesh.FindSurface(sid.c_str());
		if (surface != nullptr)
		{
			FENodeList list = surface->GetNodeList();
			for (int i = 0; i < list.Size(); ++i)
			{
				const FENode* node = list.Node(i);
				if (node) nodes.push_back(node->GetID());
			}
		}
		else
		{
			FEFacetSet* facets = mesh.FindFacetSet(sid.c_str());
			if (facets != nullptr)
			{
				FENodeList list = facets->GetNodeList();
				for (int i = 0; i < list.Size(); ++i)
				{
					const FENode* node = list.Node(i);
					if (node) nodes.push_back(node->GetID());
				}
			}
		}
	}

	std::vector<std::vector<vec3d>> forces(K, std::vector<vec3d>(T, vec3d(0, 0, 0)));
	for (size_t t = 0; t < T; ++t)
	{
		const auto& samples = loadSteps[t].loads.Samples();
		for (size_t s = 0; s < samples.size(); ++s)
		{
			for (size_t k = 0; k < K; ++k)
			{
				if (surfaces[k] == samples[s].id)
				{
					forces[k][t] = samples[s].load;
					break;
				}
			}
		}
	}

	workHistories.resize(N);
	for (size_t i = 0; i < N; ++i) workHistories[i].work.assign(T, 0.0);

	for (size_t i = 0; i < N; ++i)
	{
		const auto& field = virtualFields[i];
		for (size_t t = 0; t < T; ++t)
		{
			double w = 0.0;
			const DisplacementHistory::TimeStep* vStep = field.history.FindStepByTime(times[t], 1e-12);
			for (size_t k = 0; k < K; ++k)
			{
				const vec3d& fk = forces[k][t];
				vec3d uk(0, 0, 0);
				if (vStep != nullptr)
				{
					const std::vector<int>& nodes = surfaceNodes[k];
					for (size_t n = 0; n < nodes.size(); ++n)
					{
						std::array<double, 3> disp{};
						if (vStep->displacements.TryGet(nodes[n], disp))
						{
							uk = vec3d(disp[0], disp[1], disp[2]);
							break;
						}
					}
				}
				feLogDebugEx(m_opt.GetFEModel(), "    surface %s : force = (%g, %g, %g) u* = (%g, %g, %g)", surfaces[k].c_str(), fk.x, fk.y, fk.z, uk.x, uk.y, uk.z);
				w += fk * uk;
			}
			workHistories[i].work[t] = w;
		}
	}

	return true;
}

bool FEVFMTask::LogDiagnostics()
{
	LogSetupDiagnostics(m_opt);
	LogStressDiagnostics(m_opt);
	return true;
}

bool FEVFMTask::ExportState(const char* szfile)
{
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
	VFMExportSession session(plotPath, *m_opt.GetFEModel());
	if (!session.AddMeasuredDisplacements(m_opt.MeasuredHistory(), exportError))
	{
		feLogErrorEx(m_opt.GetFEModel(), exportError.c_str());
		return false;
	}
	if (!session.AddVirtualDisplacements(m_opt.VirtualFields(), exportError))
	{
		feLogErrorEx(m_opt.GetFEModel(), exportError.c_str());
		return false;
	}
	if (!session.AddVirtualDeformationGradients(m_opt.VirtualDeformationGradients(), exportError))
	{
		feLogErrorEx(m_opt.GetFEModel(), exportError.c_str());
		return false;
	}
	if (!session.AddMeasuredDeformationGradients(m_opt.DeformationHistory(), exportError))
	{
		feLogErrorEx(m_opt.GetFEModel(), exportError.c_str());
		return false;
	}
	if (!session.AddMeasuredStress(m_opt.StressTimeline(), exportError))
	{
		feLogErrorEx(m_opt.GetFEModel(), exportError.c_str());
		return false;
	}
	if (!session.Finalize(exportError))
	{
		feLogErrorEx(m_opt.GetFEModel(), exportError.c_str());
		return false;
	}

	feLog("VFM: exported kinematic snapshot to %s\n", plotPath.c_str());
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
