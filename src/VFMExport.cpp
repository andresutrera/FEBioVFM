#include "VFMExport.h"

#include "DisplacementContainer.h"
#include "DeformationGradientField.h"
#include "StressField.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include <FEBioPlot/FEBioPlotFile.h>
#include <FECore/writeplot.h>
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FESolidDomain.h>
#include <FECore/units.h>

namespace
{

class VFMPlotFile : public FEBioPlotFile
{
public:
    using FEBioPlotFile::FEBioPlotFile;
    using PlotFile::AddVariable;
};
// -----------------------------------------------------------------------------'

class VFMPlotMeasuredDisplacement : public FEPlotNodeData
{
public:
	VFMPlotMeasuredDisplacement(FEModel* fem)
		: FEPlotNodeData(fem, PLT_VEC3F, FMT_NODE)
	{
		SetUnits(UNIT_LENGTH);
	}

	void SetDisplacementData(const DisplacementContainer* data) { m_data = data; }

	bool Save(FEMesh& mesh, FEDataStream& a) override
	{
		const DisplacementContainer* data = m_data;
		writeNodalValues<vec3d>(mesh, a, [&](const FENode& node) {
			std::array<double, 3> disp{};
			if (data && data->TryGet(node.GetID(), disp))
			{
				return vec3d(disp[0], disp[1], disp[2]);
			}
			return vec3d(0.0, 0.0, 0.0);
		});
		return true;
	}

private:
	const DisplacementContainer* m_data = nullptr;
};

// -----------------------------------------------------------------------------'

class VFMPlotVirtualDisplacement : public FEPlotNodeData
{
public:
	VFMPlotVirtualDisplacement(FEModel* fem)
		: FEPlotNodeData(fem, PLT_VEC3F, FMT_NODE)
	{
		SetUnits(UNIT_LENGTH);
	}

	void SetDisplacementData(const DisplacementContainer* data) { m_data = data; }

	bool Save(FEMesh& mesh, FEDataStream& a) override
	{
		const DisplacementContainer* data = m_data;
		writeNodalValues<vec3d>(mesh, a, [&](const FENode& node) {
			std::array<double, 3> disp{};
			if (data && data->TryGet(node.GetID(), disp))
			{
				return vec3d(disp[0], disp[1], disp[2]);
			}
			return vec3d(0.0, 0.0, 0.0);
		});
		return true;
	}

private:
	const DisplacementContainer* m_data = nullptr;
};

// -----------------------------------------------------------------------------'

class VFMPlotDeformationGradient : public FEPlotDomainData
{
public:
	VFMPlotDeformationGradient(FEModel* fem)
		: FEPlotDomainData(fem, PLT_MAT3F, FMT_ITEM)
	{
		SetUnits(UNIT_NONE);
	}

	void SetDeformationData(const DeformationGradientField* field) { m_field = field; }

	bool Save(FEDomain& dom, FEDataStream& a) override
	{
		if (dom.Class() != FE_DOMAIN_SOLID)
		{
			for (int i = 0; i < dom.Elements(); ++i)
			{
				a << mat3d::identity();
			}
			return true;
		}

		FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
		for (int i = 0; i < sd.Elements(); ++i)
		{
			FESolidElement& el = sd.Element(i);
			const DeformationGradientField* field = m_field;
			const GaussPointDeformation* gpF = (field ? field->Find(el.GetID()) : nullptr);

			mat3d Fav = mat3d::identity();
			if (gpF && !gpF->gradients.empty())
			{
				mat3d accumulator; accumulator.zero();
				for (const mat3d& F : gpF->gradients) accumulator += F;
				accumulator /= static_cast<double>(gpF->gradients.size());
				Fav = accumulator;
			}

			a << Fav;
		}
		return true;
	}

private:
	const DeformationGradientField* m_field = nullptr;
};

// -----------------------------------------------------------------------------'

class VFMPlotStress : public FEPlotDomainData
{
public:
	VFMPlotStress(FEModel* fem)
		: FEPlotDomainData(fem, PLT_MAT3FS, FMT_ITEM)
	{
		SetUnits(UNIT_PRESSURE);
	}

	void SetStressData(const StressField* field) { m_field = field; }

	bool Save(FEDomain& dom, FEDataStream& a) override
	{
		if (dom.Class() != FE_DOMAIN_SOLID)
		{
			mat3ds zero; zero.zero();
			for (int i = 0; i < dom.Elements(); ++i) a << zero;
			return true;
		}

		FESolidDomain& sd = static_cast<FESolidDomain&>(dom);
		for (int i = 0; i < sd.Elements(); ++i)
		{
			FESolidElement& el = sd.Element(i);
			const StressField* field = m_field;
			const GaussPointStress* gpS = (field ? field->Find(el.GetID()) : nullptr);

			mat3ds Sav; Sav.zero();
			if (gpS && !gpS->stresses.empty())
			{
				mat3ds accumulator; accumulator.zero();
				for (const mat3ds& s : gpS->stresses) accumulator += s;
				accumulator /= static_cast<double>(gpS->stresses.size());
				Sav = accumulator;
			}

			a << Sav;
		}
		return true;
	}

private:
	const StressField* m_field = nullptr;
};

} // namespace

// -----------------------------------------------------------------------------'

bool ExportVFMKinematics(const std::string& filePath,
	FEModel& fem,
	const DisplacementHistory& measuredHist,
	const DisplacementHistory& virtualHist,
	const DeformationGradientHistory& defHist,
	const StressHistory& stressHist,
	std::string& error)
{
	VFMPlotFile plt(&fem);

	auto* measuredField = new VFMPlotMeasuredDisplacement(&fem);
	if (!plt.AddVariable(measuredField, "vfm_measured_disp"))
	{
		delete measuredField;
		error = "Failed to register measured displacement field.";
		return false;
	}

	auto* displacementField = new VFMPlotMeasuredDisplacement(&fem);
	if (!plt.AddVariable(displacementField, "displacement"))
	{
		delete displacementField;
		error = "Failed to register displacement field.";
		return false;
	}

	auto* virtualField = new VFMPlotVirtualDisplacement(&fem);
	if (!plt.AddVariable(virtualField, "vfm_virtual_disp"))
	{
		delete virtualField;
		error = "Failed to register virtual displacement field.";
		return false;
	}

	auto* defGradField = new VFMPlotDeformationGradient(&fem);
	if (!plt.AddVariable(defGradField, "vfm_deformation_gradient"))
	{
		delete defGradField;
		error = "Failed to register deformation gradient field.";
		return false;
	}

	auto* stressField = new VFMPlotStress(&fem);
	if (!plt.AddVariable(stressField, "vfm_stress"))
	{
		delete stressField;
		error = "Failed to register stress field.";
		return false;
	}

	plt.SetSoftwareString("FEBio VFM plug-in");

	if (!plt.Open(filePath.c_str()))
	{
		error = "Unable to create plot file: " + filePath;
		return false;
	}

	std::vector<double> times;
	times.reserve(measuredHist.Steps() + virtualHist.Steps() + defHist.Steps() + stressHist.Steps());

	const double TIME_EPS = 1e-12;

	for (const auto& step : measuredHist.StepsRef())
	{
		times.push_back(step.time);
	}

	for (const auto& step : virtualHist.StepsRef())
	{
		times.push_back(step.time);
	}

	for (const auto& step : defHist.StepsRef())
	{
		times.push_back(step.time);
	}

	for (const auto& step : stressHist.StepsRef())
	{
		times.push_back(step.time);
	}

	if (times.empty())
	{
		error = "No data available for export.";
		plt.Close();
		return false;
	}

	std::sort(times.begin(), times.end());
	std::vector<double> uniqueTimes;
	uniqueTimes.reserve(times.size());
	for (double t : times)
	{
		if (uniqueTimes.empty() || std::fabs(uniqueTimes.back() - t) > TIME_EPS)
		{
			uniqueTimes.push_back(t);
		}
	}

	auto findMeasuredStep = [&](double t) -> const DisplacementHistory::TimeStep*
	{
		return measuredHist.FindStepByTime(t, TIME_EPS);
	};

	auto findVirtualStep = [&](double t) -> const DisplacementHistory::TimeStep*
	{
		return virtualHist.FindStepByTime(t, TIME_EPS);
	};

	auto findDefStep = [&](double t) -> const DeformationGradientHistory::TimeStep*
	{
		for (const auto& step : defHist.StepsRef())
		{
			if (std::fabs(step.time - t) <= TIME_EPS) return &step;
		}
		return nullptr;
	};

	auto findStressStep = [&](double t) -> const StressHistory::TimeStep*
	{
		for (const auto& step : stressHist.StepsRef())
		{
			if (std::fabs(step.time - t) <= TIME_EPS) return &step;
		}
		return nullptr;
	};

	for (double t : uniqueTimes)
	{
		const auto* measuredStep = findMeasuredStep(t);
		const auto* virtualStep = findVirtualStep(t);
		const auto* defStep = findDefStep(t);
		const auto* stressStep = findStressStep(t);

		measuredField->SetDisplacementData(measuredStep ? &measuredStep->displacements : nullptr);
		displacementField->SetDisplacementData(measuredStep ? &measuredStep->displacements : nullptr);
		virtualField->SetDisplacementData(virtualStep ? &virtualStep->displacements : nullptr);
		defGradField->SetDeformationData(defStep ? &defStep->field : nullptr);
		stressField->SetStressData(stressStep ? &stressStep->field : nullptr);

		if (!plt.Write(static_cast<float>(t)))
		{
			error = "Failed to write plot state.";
			plt.Close();
			return false;
		}

		// feLog("VFM export: wrote plot state for t = %g\n", t);
	}

	plt.Close();
	return true;
}
