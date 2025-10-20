
#include <FEBioMech/stdafx.h>
#include "FEData.h"
#include "FEVFMInput.h"
#include "VFMStress.h"
#include "levmar/levmar.h"
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FESolidElement.h>
#include <algorithm>
#include <limits>
#include <unordered_map>

namespace
{

	double DoubleContraction(const mat3d &A, const mat3d &B)
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

	mat3d VirtualGradientFromDeformation(const mat3d &Fstar)
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
FEModelParameterVFM::FEModelParameterVFM(FEModel *fem) : FEInputParameterVFM(fem)
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
	FEModel &fem = *GetFEModel();
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
		double *pd = (double *)val.data_ptr();
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
	if (m_pd)
		return *m_pd;
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
	if (m_pd == 0)
		return false;
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
FEOptimizeDataVFM::FEOptimizeDataVFM(FEModel *fem) : m_fem(fem)
{
	// m_pSolver = 0;
	// m_pTask = 0;
	// m_niter = 0;
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
	// don't plot anything in the step
	for (int i = 0; i < m_fem->Steps(); ++i)
	{
		m_fem->GetStep(i)->SetPlotLevel(FE_PLOT_NEVER);
	}

	// initialize all input parameters
	for (int i = 0; i < (int)m_Var.size(); ++i)
	{
		FEInputParameterVFM *p = m_Var[i];
		if ((p == 0) || (p->Init() == false))
			return false;

		// set the initial value
		p->SetValue(p->InitValue());
	}

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
	if (in.Input(szfile, this) == false)
		return false;
	return true;
}

//-----------------------------------------------------------------------------
// bool FEOptimizeDataVFM::RunTask()
// {

// }

bool FEOptimizeDataVFM::SetParameterVector(const std::vector<double> &values, std::string &errorMessage)
{
	if (values.size() != m_Var.size())
	{
		errorMessage = "Parameter vector length mismatch.";
		return false;
	}

	for (size_t i = 0; i < m_Var.size(); ++i)
	{
		FEInputParameterVFM *param = m_Var[i];
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

void FEOptimizeDataVFM::GetParameterVector(std::vector<double> &values) const
{
	values.resize(m_Var.size());
	for (size_t i = 0; i < m_Var.size(); ++i)
	{
		FEInputParameterVFM *param = m_Var[i];
		values[i] = (param ? param->GetValue() : 0.0);
	}
}

bool FEOptimizeDataVFM::ComputeStress(const std::vector<double> &params, std::string &errorMessage)
{

	if (!SetParameterVector(params, errorMessage))
		return false;

	auto &defHistory = DeformationHistory();
	auto &stressHistory = StressTimeline();
	auto &piolaHistory = FirstPiolaTimeline();

	stressHistory.Clear();
	stressHistory.Reserve(defHistory.Steps());
	piolaHistory.Clear();
	piolaHistory.Reserve(defHistory.Steps());

	if (defHistory.Empty())
		return true;

	for (const auto &defStep : defHistory)
	{
		auto &stressStep = stressHistory.AddStep(defStep.time);
		auto &piolaStep = piolaHistory.AddStep(defStep.time);

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

std::vector<double> FEOptimizeDataVFM::ComputeInternalWork(const std::vector<double> &params)
{
	// set parameters
	std::string _;
	if (!SetParameterVector(params, _))
		return {};

	// ensure stress/first-Piola are current
	if (!ComputeStress(params, _))
		return {};

	const size_t fieldCount = m_virtualDefGradients.Size();
	const size_t timeCount = FirstPiolaTimeline().Steps();
	if (fieldCount == 0 || timeCount == 0)
		return {};

	// precompute integration weights per element
	FEMesh &mesh = GetFEModel()->GetMesh();
	std::unordered_map<int, std::vector<double>> wmap;
	wmap.reserve(static_cast<size_t>(mesh.Domains()) * 8);

	for (int d = 0; d < mesh.Domains(); ++d)
	{
		auto *sd = dynamic_cast<FESolidDomain *>(&mesh.Domain(d));
		if (!sd)
			continue;
		for (int e = 0; e < sd->Elements(); ++e)
		{
			auto &el = static_cast<FESolidElement &>(sd->ElementRef(e));
			const int nint = el.GaussPoints();
			std::vector<double> w(nint, 0.0);
			double *gw = el.GaussWeights();
			for (int n = 0; n < nint; ++n)
			{
				const double gwgt = gw ? gw[n] : 1.0;
				w[n] = gwgt * sd->detJ0(el, n);
			}
			wmap.emplace(el.GetID(), std::move(w));
		}
	}

	std::vector<double> internalWork(fieldCount * timeCount, 0.0);
	const double tTol = 1e-12;

	for (size_t f = 0; f < fieldCount; ++f)
	{
		const auto &vfield = m_virtualDefGradients[f];
		for (size_t t = 0; t < timeCount; ++t)
		{
			const auto &Pstep = m_firstPiolaHistory[t];
			const auto *Vstep = vfield.history.FindStepByTime(Pstep.time, tTol);

			double iw = 0.0;
			for (const auto &gpP : Pstep.field.Data())
			{
				auto it = wmap.find(gpP.elementId);
				if (it == wmap.end())
					return {}; // silent fail per request
				const auto &w = it->second;

				const auto *gpV = Vstep ? Vstep->field.Find(gpP.elementId) : nullptr;
				const size_t nint = gpP.stresses.size();

				for (size_t n = 0; n < nint; ++n)
				{
					const mat3d &P = gpP.stresses[n];
					mat3d G;
					G.zero();
					if (gpV && n < gpV->gradients.size())
					{
						G = VirtualGradientFromDeformation(gpV->gradients[n]);
					}
					iw += DoubleContraction(P, G) * w[n];
				}
			}
			internalWork[f * timeCount + t] = iw;
		}
	}
	return internalWork;
}