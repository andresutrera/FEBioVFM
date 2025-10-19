
#include <FEBioMech/stdafx.h>
#include "FEData.h"
#include "FEVFMInput.h"
#include "VFMStress.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FESolidElement.h>
#include <unordered_map>

namespace
{

double DoubleContraction(const mat3d& A, const mat3d& B)
{
	double s = 0.0;
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			s += A[i][j] * B[i][j];
		}
	}
	return s;
}

mat3d VirtualGradientFromDeformation(const mat3d& Fstar)
{
	mat3d G = Fstar;
	G[0][0] -= 1.0;
	G[1][1] -= 1.0;
	G[2][2] -= 1.0;
	return G;
}

} // namespace
/**
 * @brief Construct a model parameter adapter with a null data pointer.
 *
 * @param fem FEBio model provided by the plugin.
 * @note The actual parameter lookup occurs inside Init() to ensure that the
 * model has fully parsed its input file before we attempt to resolve names.
 */
FEModelParameterVFM::FEModelParameterVFM(FEModel* fem) : FEInputParameterVFM(fem)
{
	m_pd = 0;
}

/**
 * @brief Resolve the FEModel parameter referenced by this adapter.
 *
 * The function uses FEModel::GetParameterValue with the stored name, verifies
 * that the parameter is a scalar double, and caches a pointer to the raw data.
 *
 * @return true if the parameter was found and is writable.
 * @note Detailed log messages explain the failure reason when the lookup fails,
 * which helps end users diagnose typos in the XML configuration.
 */
bool FEModelParameterVFM::Init()
{
	// find the variable
	FEModel& fem = *GetFEModel();
	string name = GetName();
	FEParamValue val = fem.GetParameterValue(ParamString(name.c_str()));

	// see if we found the parameter
	if (val.isValid() == false)
	{
		feLogError("Cannot find parameter %s", name.c_str());
		return false;
	}

	// Make sure it's a double
	if (val.type() == FE_PARAM_DOUBLE)
	{
		// make sure we have a valid data pointer
		double* pd = (double*)val.data_ptr();
		if (pd == 0)
		{
			feLogError("Invalid data pointer for parameter %s", name.c_str());
			return false;
		}

		// store the pointer to the parameter
		m_pd = pd;
	}
	else
	{
		feLogError("Invalid parameter type for parameter %s", name.c_str());
		return false;
	}

	return true;
}

/**
 * @brief Get the current value of the FEBio model parameter.
 *
 * @return The parameter value when the pointer is valid, otherwise zero.
 * @note Returning zero on failure matches the legacy plugin behaviour, keeping
 * the scaffold compatible with existing scripts that assume this fallback.
 */
double FEModelParameterVFM::GetValue()
{
	if (m_pd) return *m_pd;
	return 0;
}

/**
 * @brief Update the FEBio model parameter with a new value.
 *
 * @param newValue Value supplied by the optimizer.
 * @return true if the cached pointer is valid and the assignment succeeded.
 * @note The method does not enforce bounds; callers are expected to rely on the
 * metadata stored in FEInputParameterVFM to vet proposed updates.
 */
bool FEModelParameterVFM::SetValue(double newValue)
{
	if (m_pd == 0) return false;
	(*m_pd) = newValue;
	return true;
}

//=============================================================================

/**
 * @brief Construct the optimization data container.
 *
 * @param fem FEBio model that participates in the Virtual Fields Method.
 * @note Several member pointers remain commented out to document the intended
 * design (solver, task, objective function) without introducing unused members.
 */
FEOptimizeDataVFM::FEOptimizeDataVFM(FEModel* fem) : m_fem(fem)
{
	// m_pSolver = 0;
	// m_pTask = 0;
	m_niter = 0;
	// m_obj = 0;
}

/**
 * @brief Destroy the optimization data container.
 *
 * @note Solver ownership is deferred; the explicit delete remains commented
 * out because the corresponding pointer is not allocated yet.
 */
FEOptimizeDataVFM::~FEOptimizeDataVFM(void)
{
	// delete m_pSolver;
}

/**
 * @brief Initialize optimization-related data prior to running the task.
 *
 * The method disables FEBio plot output for all load steps, blocks logging while
 * placeholder task initialization would occur, and asks each registered input
 * parameter to initialize itself and apply its starting value.
 *
 * @return true if all registered parameters initialized successfully.
 * @note Objective function and solver initialization hooks are still TODO and
 * appear as commented code to preserve the roadmap in the source and docs.
 */
bool FEOptimizeDataVFM::Init()
{
	// allocate default optimization solver if none specified in input file
	// if (m_pSolver == 0) m_pSolver = new FELMOptimizeMethod(GetFEModel());

	// allocate default solver if none specified in input file
	// if (m_pTask == 0) m_pTask = fecore_new<FECoreTask>("solve", m_fem);

	// don't plot anything
	for (int i = 0; i < m_fem->Steps(); ++i)
	{
		m_fem->GetStep(i)->SetPlotLevel(FE_PLOT_NEVER);
	}

	// do the initialization of the task
	GetFEModel()->BlockLog();
	// if (m_pTask->Init(0) == false) return false;
	GetFEModel()->UnBlockLog();

	m_initialParameters.clear();
	m_initialParameters.resize(m_Var.size(), 0.0);

	// initialize all input parameters
	for (int i=0; i<(int)m_Var.size(); ++i)
	{
		FEInputParameterVFM* p = m_Var[i];
		if ((p==0) || (p->Init() == false)) return false;

		// set the initial value
		p->SetValue(p->InitValue());
		m_initialParameters[i] = p->GetValue();
	}

	// initialize the objective function
	// if (m_obj == 0) return false;
	// if (m_obj->Init() == false) return false;

	return true;
}

/**
 * @brief Placeholder for the optimization loop.
 *
 * @return false until the routine is implemented.
 * @note The previous implementation was not available; returning false makes
 * the incomplete status explicit to both callers and the generated documentation.
 */
bool FEOptimizeDataVFM::Solve()
{
	return false;
}

/**
 * @brief Read and process the VFM-specific optimization input file.
 *
 * @param szfile Path to the XML input file.
 * @return true when parsing succeeded.
 */
bool FEOptimizeDataVFM::Input(const char *szfile)
{
	FEVFMInput in;
	if (in.Input(szfile, this) == false) return false;
	return true;
}

//-----------------------------------------------------------------------------
// bool FEOptimizeDataVFM::RunTask()
// {

// }

/**
 * @brief Solve the FE problem with a new set of optimization parameters.
 *
 * @param a Vector containing trial parameter values.
 * @return false until the implementation is provided.
 * @note The method is intentionally left as a stub because the forward solve is
 * tightly coupled to FEBio infrastructure that has not been integrated yet.
 */
bool FEOptimizeDataVFM::FESolve(const vector<double>& a)
{
	return false;
}

bool FEOptimizeDataVFM::SetParameterVector(const std::vector<double>& values, std::string& errorMessage)
{
	if (values.size() != m_Var.size())
	{
		errorMessage = "Parameter vector length mismatch.";
		return false;
	}

	for (size_t i = 0; i < m_Var.size(); ++i)
	{
		FEInputParameterVFM* param = m_Var[i];
		if (param == nullptr)
		{
			errorMessage = "Encountered an uninitialized parameter slot at index " + std::to_string(i) + ".";
			return false;
		}

		if (!param->SetValue(values[i]))
		{
			errorMessage = "Failed to assign parameter '" + param->GetName() + "'.";
			return false;
		}
	}

	return true;
}

void FEOptimizeDataVFM::GetParameterVector(std::vector<double>& values) const
{
	values.resize(m_Var.size());
	for (size_t i = 0; i < m_Var.size(); ++i)
	{
		FEInputParameterVFM* param = m_Var[i];
		values[i] = (param ? param->GetValue() : 0.0);
	}
}

bool FEOptimizeDataVFM::ResetParametersToInitial(std::string& errorMessage)
{
	if (m_initialParameters.size() != m_Var.size())
	{
		errorMessage = "Initial parameter vector was not captured correctly.";
		return false;
	}
	return SetParameterVector(m_initialParameters, errorMessage);
}

bool FEOptimizeDataVFM::RebuildStressHistories(std::string& errorMessage)
{
	return RebuildStressHistoriesInternal(errorMessage);
}

bool FEOptimizeDataVFM::RebuildStressHistories(const std::vector<double>& parameterValues,
	bool restoreOriginalValues,
	std::string& errorMessage)
{
	std::vector<double> originalValues;
	if (restoreOriginalValues)
	{
		GetParameterVector(originalValues);
	}

	if (!parameterValues.empty())
	{
		if (!SetParameterVector(parameterValues, errorMessage))
		{
			return false;
		}
	}

	const bool result = RebuildStressHistoriesInternal(errorMessage);

	if (restoreOriginalValues)
	{
		std::string restoreError;
		if (!SetParameterVector(originalValues, restoreError))
		{
			errorMessage = restoreError;
			return false;
		}
		if (!RebuildStressHistoriesInternal(restoreError))
		{
			errorMessage = restoreError;
			return false;
		}
	}

	return result;
}

bool FEOptimizeDataVFM::RebuildStressHistoriesInternal(std::string& errorMessage)
{
	auto& defHistory = DeformationHistory();
	auto& stressHistory = StressTimeline();
	auto& piolaHistory = FirstPiolaTimeline();

	stressHistory.Clear();
	stressHistory.Reserve(defHistory.Steps());
	piolaHistory.Clear();
	piolaHistory.Reserve(defHistory.Steps());

	if (defHistory.Empty()) return true;

	for (const auto& defStep : defHistory)
	{
		auto& stressStep = stressHistory.AddStep(defStep.time);
		auto& piolaStep = piolaHistory.AddStep(defStep.time);

		if (!VFMStress::ComputeCauchyStress(*GetFEModel(),
			defStep.field,
			stressStep.field,
			errorMessage))
		{
			return false;
		}

		if (!VFMStress::ComputeFirstPiolaStress(defStep.field,
			stressStep.field,
			piolaStep.field,
			errorMessage))
		{
			return false;
		}
	}

	return true;
}

bool FEOptimizeDataVFM::AssembleResidual(std::vector<double>& residual, std::string& errorMessage)
{
	return AssembleResidualInternal(residual, errorMessage);
}

bool FEOptimizeDataVFM::AssembleResidual(const std::vector<double>& parameterValues,
	bool restoreOriginalValues,
	std::vector<double>& residual,
	std::string& errorMessage)
{
	std::vector<double> originalValues;
	if (restoreOriginalValues)
	{
		GetParameterVector(originalValues);
	}

	bool parametersUpdated = false;
	if (!parameterValues.empty())
	{
		if (!SetParameterVector(parameterValues, errorMessage))
		{
			return false;
		}
		parametersUpdated = true;
	}

	if (parametersUpdated)
	{
		if (!RebuildStressHistoriesInternal(errorMessage))
		{
			if (restoreOriginalValues)
			{
				std::string restoreError;
				if (SetParameterVector(originalValues, restoreError))
				{
					if (!RebuildStressHistoriesInternal(restoreError))
					{
						errorMessage = restoreError;
					}
				}
				else
				{
					errorMessage = restoreError;
				}
			}
			return false;
		}
	}

	const bool assembled = AssembleResidualInternal(residual, errorMessage);

	if (restoreOriginalValues)
	{
		std::string restoreError;
		if (!SetParameterVector(originalValues, restoreError))
		{
			errorMessage = restoreError;
			return false;
		}
		if (!RebuildStressHistoriesInternal(restoreError))
		{
			errorMessage = restoreError;
			return false;
		}
	}

	return assembled;
}

bool FEOptimizeDataVFM::AssembleResidualInternal(std::vector<double>& residual, std::string& errorMessage)
{
	const size_t fieldCount = m_virtualDefGradients.Size();
	const size_t timeCount = FirstPiolaTimeline().Steps();

	residual.assign(fieldCount * timeCount, 0.0);

	if ((fieldCount == 0) || (timeCount == 0))
	{
		residual.clear();
		return true;
	}

	if (m_virtualExternalWork.size() != fieldCount)
	{
		errorMessage = "Virtual external work history count does not match the number of virtual fields.";
		return false;
	}

	for (size_t i = 0; i < fieldCount; ++i)
	{
		if (m_virtualExternalWork[i].work.size() < timeCount)
		{
			errorMessage = "Virtual external work history for field index " + std::to_string(i) + " is shorter than the stress timeline.";
			return false;
		}
	}

	FEMesh& mesh = GetFEModel()->GetMesh();
	std::unordered_map<int, std::vector<double>> integrationWeights;

	const int domainCount = mesh.Domains();
	for (int domIdx = 0; domIdx < domainCount; ++domIdx)
	{
		FEDomain& domain = mesh.Domain(domIdx);
		auto* solidDomain = dynamic_cast<FESolidDomain*>(&domain);
		if (solidDomain == nullptr) continue;

		const int elemCount = solidDomain->Elements();
		for (int elemIdx = 0; elemIdx < elemCount; ++elemIdx)
		{
			FESolidElement& el = static_cast<FESolidElement&>(solidDomain->ElementRef(elemIdx));
			const int elemId = el.GetID();
			const int nint = el.GaussPoints();

			std::vector<double> weights(nint, 0.0);
			double* gw = el.GaussWeights();
			for (int n = 0; n < nint; ++n)
			{
				const double gaussWeight = (gw ? gw[n] : 1.0);
				const double detJ0 = solidDomain->detJ0(el, n);
				weights[n] = gaussWeight * detJ0;
			}

			integrationWeights.emplace(elemId, std::move(weights));
		}
	}

	const double timeTolerance = 1e-12;
	for (size_t fieldIdx = 0; fieldIdx < fieldCount; ++fieldIdx)
	{
		const auto& virtualField = m_virtualDefGradients[fieldIdx];
		const auto& externalWork = m_virtualExternalWork[fieldIdx].work;

		for (size_t timeIdx = 0; timeIdx < timeCount; ++timeIdx)
		{
			const FirstPiolaHistory::TimeStep& piolaStep = m_firstPiolaHistory[timeIdx];
			const DeformationGradientHistory::TimeStep* virtualStep = virtualField.history.FindStepByTime(piolaStep.time, timeTolerance);

			double internalVirtualWork = 0.0;
			for (const GaussPointFirstPiola& gpPiola : piolaStep.field.Data())
			{
				auto weightIt = integrationWeights.find(gpPiola.elementId);
				if (weightIt == integrationWeights.end())
				{
					errorMessage = "Missing integration weights for element " + std::to_string(gpPiola.elementId) + ".";
					return false;
				}

				const std::vector<double>& weights = weightIt->second;
				const GaussPointDeformation* gpVirtual = (virtualStep ? virtualStep->field.Find(gpPiola.elementId) : nullptr);

				const size_t gaussCount = gpPiola.stresses.size();
				if (weights.size() < gaussCount)
				{
					errorMessage = "Integration weight count mismatch for element " + std::to_string(gpPiola.elementId) + ".";
					return false;
				}

				for (size_t n = 0; n < gaussCount; ++n)
				{
					const mat3d& P = gpPiola.stresses[n];
					mat3d G; G.zero();
					if (gpVirtual && n < gpVirtual->gradients.size())
					{
						G = VirtualGradientFromDeformation(gpVirtual->gradients[n]);
					}

					internalVirtualWork += DoubleContraction(P, G) * weights[n];
				}
			}

			const size_t residualIndex = fieldIdx * timeCount + timeIdx;
			residual[residualIndex] = internalVirtualWork - externalWork[timeIdx];
		}
	}

	return true;
}
