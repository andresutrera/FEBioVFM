#include "VFMExport.h"

#include "DisplacementContainer.h"
#include "DeformationGradientField.h"
#include "StressField.h"
#include "FirstPiolaField.h"
#include "FEData.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include <FEBioPlot/FEBioPlotFile.h>
#include <FECore/writeplot.h>
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FESolidDomain.h>
#include <FECore/units.h>
#include <FECore/log.h>

namespace
{

	class VFMPlotFile : public FEBioPlotFile
	{
	public:
		using FEBioPlotFile::FEBioPlotFile;
		using PlotFile::AddVariable;
	};

	class VFMPlotMeasuredDisplacement : public FEPlotNodeData
	{
	public:
		explicit VFMPlotMeasuredDisplacement(FEModel *fem)
			: FEPlotNodeData(fem, PLT_VEC3F, FMT_NODE)
		{
			SetUnits(UNIT_LENGTH);
		}

		void SetDisplacementData(const DisplacementContainer *data) { m_data = data; }

		bool Save(FEMesh &mesh, FEDataStream &a) override
		{
			const DisplacementContainer *data = m_data;
			writeNodalValues<vec3d>(mesh, a, [&](const FENode &node)
									{
			std::array<double, 3> disp{};
			if (data && data->TryGet(node.GetID(), disp))
			{
				return vec3d(disp[0], disp[1], disp[2]);
			}
			return vec3d(0.0, 0.0, 0.0); });
			return true;
		}

	private:
		const DisplacementContainer *m_data = nullptr;
	};

	class VFMPlotVirtualDisplacement : public FEPlotNodeData
	{
	public:
		explicit VFMPlotVirtualDisplacement(FEModel *fem)
			: FEPlotNodeData(fem, PLT_VEC3F, FMT_NODE)
		{
			SetUnits(UNIT_LENGTH);
		}

		void SetDisplacementData(const DisplacementContainer *data) { m_data = data; }

		bool Save(FEMesh &mesh, FEDataStream &a) override
		{
			const DisplacementContainer *data = m_data;
			writeNodalValues<vec3d>(mesh, a, [&](const FENode &node)
									{
			std::array<double, 3> disp{};
			if (data && data->TryGet(node.GetID(), disp))
			{
				return vec3d(disp[0], disp[1], disp[2]);
			}
			return vec3d(0.0, 0.0, 0.0); });
			return true;
		}

	private:
		const DisplacementContainer *m_data = nullptr;
	};

	class VFMPlotDeformationGradient : public FEPlotDomainData
	{
	public:
		explicit VFMPlotDeformationGradient(FEModel *fem)
			: FEPlotDomainData(fem, PLT_MAT3F, FMT_ITEM)
		{
			SetUnits(UNIT_NONE);
		}

		void SetDeformationData(const DeformationGradientField *field) { m_field = field; }

		bool Save(FEDomain &dom, FEDataStream &a) override
		{
			if (dom.Class() != FE_DOMAIN_SOLID)
			{
				for (int i = 0; i < dom.Elements(); ++i)
				{
					a << mat3d::identity();
				}
				return true;
			}

			FESolidDomain &sd = static_cast<FESolidDomain &>(dom);
			for (int i = 0; i < sd.Elements(); ++i)
			{
				FESolidElement &el = sd.Element(i);
				const DeformationGradientField *field = m_field;
				const GaussPointDeformation *gpF = (field ? field->Find(el.GetID()) : nullptr);

				mat3d Fav = mat3d::identity();
				if (gpF && !gpF->gradients.empty())
				{
					mat3d accumulator;
					accumulator.zero();
					for (const mat3d &F : gpF->gradients)
						accumulator += F;
					accumulator /= static_cast<double>(gpF->gradients.size());
					Fav = accumulator;
				}

				a << Fav;
			}
			return true;
		}

	private:
		const DeformationGradientField *m_field = nullptr;
	};

	class VFMPlotStress : public FEPlotDomainData
	{
	public:
		explicit VFMPlotStress(FEModel *fem)
			: FEPlotDomainData(fem, PLT_MAT3FS, FMT_ITEM)
		{
			SetUnits(UNIT_PRESSURE);
		}

		void SetStressData(const StressField *field) { m_field = field; }

		bool Save(FEDomain &dom, FEDataStream &a) override
		{
			if (dom.Class() != FE_DOMAIN_SOLID)
			{
				mat3ds zero;
				zero.zero();
				for (int i = 0; i < dom.Elements(); ++i)
					a << zero;
				return true;
			}

			FESolidDomain &sd = static_cast<FESolidDomain &>(dom);
			for (int i = 0; i < sd.Elements(); ++i)
			{
				FESolidElement &el = sd.Element(i);
				const StressField *field = m_field;
				const GaussPointStress *gpS = (field ? field->Find(el.GetID()) : nullptr);

				mat3ds Sav;
				Sav.zero();
				if (gpS && !gpS->stresses.empty())
				{
					mat3ds accumulator;
					accumulator.zero();
					for (const mat3ds &s : gpS->stresses)
						accumulator += s;
					accumulator /= static_cast<double>(gpS->stresses.size());
					Sav = accumulator;
				}

				a << Sav;
			}
			return true;
		}

	private:
		const StressField *m_field = nullptr;
	};

	class VFMPlotFirstPiolaStress : public FEPlotDomainData
	{
	public:
		explicit VFMPlotFirstPiolaStress(FEModel *fem)
			: FEPlotDomainData(fem, PLT_MAT3F, FMT_ITEM)
		{
			SetUnits(UNIT_PRESSURE);
		}

		void SetStressData(const FirstPiolaField *field) { m_field = field; }

		bool Save(FEDomain &dom, FEDataStream &a) override
		{
			if (dom.Class() != FE_DOMAIN_SOLID)
			{
				mat3d zero;
				zero.zero();
				for (int i = 0; i < dom.Elements(); ++i)
					a << zero;
				return true;
			}

			FESolidDomain &sd = static_cast<FESolidDomain &>(dom);
			for (int i = 0; i < sd.Elements(); ++i)
			{
				FESolidElement &el = sd.Element(i);
				const FirstPiolaField *field = m_field;
				const GaussPointFirstPiola *gpP = (field ? field->Find(el.GetID()) : nullptr);

				mat3d Pav;
				Pav.zero();
				if (gpP && !gpP->stresses.empty())
				{
					mat3d accumulator;
					accumulator.zero();
					for (const mat3d &P : gpP->stresses)
						accumulator += P;
					accumulator /= static_cast<double>(gpP->stresses.size());
					Pav = accumulator;
				}

				a << Pav;
			}
			return true;
		}

	private:
		const FirstPiolaField *m_field = nullptr;
	};

} // namespace

struct VFMExportSession::Impl
{
	FEModel &fem;
	std::string filePath;
	VFMPlotFile plot;

	VFMPlotMeasuredDisplacement *measuredField = nullptr;
	VFMPlotMeasuredDisplacement *measuredEvalField = nullptr;
	const DisplacementHistory *measuredHist = nullptr;

	std::vector<VFMPlotVirtualDisplacement *> virtualPlots;
	std::vector<const VirtualDisplacementCollection::Field *> virtualRefs;

	std::vector<VFMPlotDeformationGradient *> virtualDefPlots;
	std::vector<const VirtualDeformationGradientCollection::Field *> virtualDefRefs;

	VFMPlotDeformationGradient *measuredDefField = nullptr;
	const DeformationGradientHistory *measuredDefHist = nullptr;

	VFMPlotStress *stressField = nullptr;
	const StressHistory *stressHist = nullptr;
	VFMPlotFirstPiolaStress *firstPiolaField = nullptr;
	const FirstPiolaHistory *firstPiolaHist = nullptr;

	std::vector<double> times;
	bool finalized = false;

	static constexpr double TIME_EPS = 1e-12;

	Impl(const std::string &path, FEModel &femRef)
		: fem(femRef), filePath(path), plot(&femRef)
	{
		plot.SetSoftwareString("FEBio VFM plug-in");
	}

	template <typename History>
	void AppendTimes(const History &history)
	{
		for (const auto &step : history.StepsRef())
		{
			times.push_back(step.time);
		}
	}

	bool AddMeasuredDisplacements(const DisplacementHistory &hist, std::string &error)
	{
		if (measuredHist != nullptr)
		{
			error = "Measured displacements already registered.";
			return false;
		}

		VFMPlotMeasuredDisplacement *measured = new VFMPlotMeasuredDisplacement(&fem);
		if (!plot.AddVariable(measured, "Measured Displacement"))
		{
			error = "Failed to register measured displacement field.";
			delete measured;
			return false;
		}

		VFMPlotMeasuredDisplacement *eval = new VFMPlotMeasuredDisplacement(&fem);
		if (!plot.AddVariable(eval, "displacement"))
		{
			error = "Failed to register displacement field.";
			delete eval;
			return false;
		}

		measuredField = measured;
		measuredEvalField = eval;
		measuredHist = &hist;
		AppendTimes(hist);
		return true;
	}

	bool AddVirtualDisplacements(const VirtualDisplacementCollection &fields, std::string &error)
	{
		if (fields.Empty())
		{
			VFMPlotVirtualDisplacement *plotVar = new VFMPlotVirtualDisplacement(&fem);
			if (!plot.AddVariable(plotVar, "Virtual Displacement"))
			{
				error = "Failed to register virtual displacement field.";
				delete plotVar;
				return false;
			}
			virtualRefs.push_back(nullptr);
			virtualPlots.push_back(plotVar);
			return true;
		}

		virtualPlots.reserve(virtualPlots.size() + fields.Size());
		virtualRefs.reserve(virtualRefs.size() + fields.Size());

		size_t idx = 0;
		for (const auto &field : fields.Data())
		{
			VFMPlotVirtualDisplacement *plotVar = new VFMPlotVirtualDisplacement(&fem);
			std::string varName = "Virtual Displacement";
			if (!field.id.empty())
			{
				varName += " " + field.id;
			}
			else if (fields.Size() > 1)
			{
				varName += " #" + std::to_string(idx);
			}

			if (!plot.AddVariable(plotVar, varName.c_str()))
			{
				error = "Failed to register virtual displacement field.";
				delete plotVar;
				return false;
			}

			virtualRefs.push_back(&field);
			virtualPlots.push_back(plotVar);
			AppendTimes(field.history);
			++idx;
		}
		return true;
	}

	bool AddVirtualDeformationGradients(const VirtualDeformationGradientCollection &fields, std::string &error)
	{
		if (fields.Empty())
			return true;

		virtualDefPlots.reserve(virtualDefPlots.size() + fields.Size());
		virtualDefRefs.reserve(virtualDefRefs.size() + fields.Size());

		size_t idx = 0;
		for (const auto &field : fields.Data())
		{
			VFMPlotDeformationGradient *plotVar = new VFMPlotDeformationGradient(&fem);
			std::string varName = "Virtual Defgrad";
			if (!field.id.empty())
			{
				varName += " " + field.id;
			}
			else if (fields.Size() > 1)
			{
				varName += " #" + std::to_string(idx);
			}

			if (!plot.AddVariable(plotVar, varName.c_str()))
			{
				error = "Failed to register virtual deformation gradient field.";
				delete plotVar;
				return false;
			}

			virtualDefRefs.push_back(&field);
			virtualDefPlots.push_back(plotVar);
			AppendTimes(field.history);
			++idx;
		}
		return true;
	}

	bool AddMeasuredDeformationGradients(const DeformationGradientHistory &hist, std::string &error)
	{
		if (measuredDefHist != nullptr)
		{
			error = "Measured deformation gradients already registered.";
			return false;
		}

		VFMPlotDeformationGradient *plotVar = new VFMPlotDeformationGradient(&fem);
		if (!plot.AddVariable(plotVar, "Measured Deformation Gradient"))
		{
			error = "Failed to register deformation gradient field.";
			delete plotVar;
			return false;
		}

		measuredDefField = plotVar;
		measuredDefHist = &hist;
		AppendTimes(hist);
		return true;
	}

	bool AddMeasuredStress(const StressHistory &hist, std::string &error)
	{
		if (stressHist != nullptr)
		{
			error = "Measured stress history already registered.";
			return false;
		}

		VFMPlotStress *plotVar = new VFMPlotStress(&fem);
		if (!plot.AddVariable(plotVar, "Estimated Stress"))
		{
			error = "Failed to register stress field.";
			delete plotVar;
			return false;
		}

		stressField = plotVar;
		stressHist = &hist;
		AppendTimes(hist);
		return true;
	}

	bool AddFirstPiolaStress(const FirstPiolaHistory &hist, std::string &error)
	{
		if (firstPiolaHist != nullptr)
		{
			error = "Measured first Piola stress history already registered.";
			return false;
		}

		VFMPlotFirstPiolaStress *plotVar = new VFMPlotFirstPiolaStress(&fem);
		if (!plot.AddVariable(plotVar, "Estimated First Piola Stress"))
		{
			error = "Failed to register first Piola stress field.";
			delete plotVar;
			return false;
		}

		firstPiolaField = plotVar;
		firstPiolaHist = &hist;
		AppendTimes(hist);
		return true;
	}

	bool Finalize(std::string &error)
	{
		if (finalized)
			return true;

		if (times.empty())
		{
			error = "No data available for export.";
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

		if (!plot.Open(filePath.c_str()))
		{
			error = "Unable to create plot file: " + filePath;
			return false;
		}

		for (double t : uniqueTimes)
		{
			const DisplacementHistory::TimeStep *measuredStep = nullptr;
			if (measuredHist)
			{
				measuredStep = measuredHist->FindStepByTime(t, TIME_EPS);
			}

			if (measuredField)
			{
				measuredField->SetDisplacementData(measuredStep ? &measuredStep->displacements : nullptr);
			}
			if (measuredEvalField)
			{
				measuredEvalField->SetDisplacementData(measuredStep ? &measuredStep->displacements : nullptr);
			}

			for (size_t i = 0; i < virtualPlots.size(); ++i)
			{
				const DisplacementHistory::TimeStep *stepPtr = nullptr;
				if (virtualRefs[i] != nullptr)
				{
					stepPtr = virtualRefs[i]->history.FindStepByTime(t, TIME_EPS);
				}
				virtualPlots[i]->SetDisplacementData(stepPtr ? &stepPtr->displacements : nullptr);
			}

			const DeformationGradientHistory::TimeStep *measuredDefStep = nullptr;
			if (measuredDefHist)
			{
				for (const auto &step : measuredDefHist->StepsRef())
				{
					if (std::fabs(step.time - t) <= TIME_EPS)
					{
						measuredDefStep = &step;
						break;
					}
				}
			}

			for (size_t i = 0; i < virtualDefPlots.size(); ++i)
			{
				const DeformationGradientHistory::TimeStep *stepPtr = nullptr;
				if (virtualDefRefs[i] != nullptr)
				{
					for (const auto &step : virtualDefRefs[i]->history.StepsRef())
					{
						if (std::fabs(step.time - t) <= TIME_EPS)
						{
							stepPtr = &step;
							break;
						}
					}
				}
				virtualDefPlots[i]->SetDeformationData(stepPtr ? &stepPtr->field : nullptr);
			}

			if (measuredDefField)
			{
				measuredDefField->SetDeformationData(measuredDefStep ? &measuredDefStep->field : nullptr);
			}

			const StressHistory::TimeStep *stressStep = nullptr;
			if (stressHist)
			{
				for (const auto &step : stressHist->StepsRef())
				{
					if (std::fabs(step.time - t) <= TIME_EPS)
					{
						stressStep = &step;
						break;
					}
				}
			}

			if (stressField)
			{
				stressField->SetStressData(stressStep ? &stressStep->field : nullptr);
			}

			const FirstPiolaHistory::TimeStep *piolaStep = nullptr;
			if (firstPiolaHist)
			{
				for (const auto &step : firstPiolaHist->StepsRef())
				{
					if (std::fabs(step.time - t) <= TIME_EPS)
					{
						piolaStep = &step;
						break;
					}
				}
			}

			if (firstPiolaField)
			{
				firstPiolaField->SetStressData(piolaStep ? &piolaStep->field : nullptr);
			}

			if (!plot.Write(static_cast<float>(t)))
			{
				error = "Failed to write plot state.";
				plot.Close();
				return false;
			}
		}

		plot.Close();
		finalized = true;
		return true;
	}
};

VFMExportSession::VFMExportSession(const std::string &filePath, FEModel &fem)
	: m_impl(std::unique_ptr<Impl>(new Impl(filePath, fem)))
{
}

VFMExportSession::~VFMExportSession() = default;

bool VFMExportSession::AddMeasuredDisplacements(const DisplacementHistory &hist, std::string &error)
{
	return m_impl->AddMeasuredDisplacements(hist, error);
}

bool VFMExportSession::AddVirtualDisplacements(const VirtualDisplacementCollection &fields, std::string &error)
{
	return m_impl->AddVirtualDisplacements(fields, error);
}

bool VFMExportSession::AddVirtualDeformationGradients(const VirtualDeformationGradientCollection &fields, std::string &error)
{
	return m_impl->AddVirtualDeformationGradients(fields, error);
}

bool VFMExportSession::AddMeasuredDeformationGradients(const DeformationGradientHistory &hist, std::string &error)
{
	return m_impl->AddMeasuredDeformationGradients(hist, error);
}

bool VFMExportSession::AddMeasuredStress(const StressHistory &hist, std::string &error)
{
	return m_impl->AddMeasuredStress(hist, error);
}

bool VFMExportSession::AddFirstPiolaStress(const FirstPiolaHistory &hist, std::string &error)
{
	return m_impl->AddFirstPiolaStress(hist, error);
}

bool VFMExportSession::Finalize(std::string &error)
{
	return m_impl->Finalize(error);
}
