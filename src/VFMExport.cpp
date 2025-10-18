#include "VFMExport.h"

#include "DisplacementContainer.h"
#include "DeformationGradientField.h"

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
	VFMPlotMeasuredDisplacement(FEModel* fem, const DisplacementContainer& data)
		: FEPlotNodeData(fem, PLT_VEC3F, FMT_NODE), m_data(data)
	{
		SetUnits(UNIT_LENGTH);
	}

	bool Save(FEMesh& mesh, FEDataStream& a) override
	{
		writeNodalValues<vec3d>(mesh, a, [&](const FENode& node) {
			std::array<double, 3> disp{};
			if (m_data.TryGet(node.GetID(), disp))
			{
				return vec3d(disp[0], disp[1], disp[2]);
			}
			return vec3d(0.0, 0.0, 0.0);
		});
		return true;
	}

private:
	const DisplacementContainer& m_data;
};

// -----------------------------------------------------------------------------'

class VFMPlotVirtualDisplacement : public FEPlotNodeData
{
public:
	VFMPlotVirtualDisplacement(FEModel* fem, const DisplacementContainer& data)
		: FEPlotNodeData(fem, PLT_VEC3F, FMT_NODE), m_data(data)
	{
		SetUnits(UNIT_LENGTH);
	}

	bool Save(FEMesh& mesh, FEDataStream& a) override
	{
		writeNodalValues<vec3d>(mesh, a, [&](const FENode& node) {
			std::array<double, 3> disp{};
			if (m_data.TryGet(node.GetID(), disp))
			{
				return vec3d(disp[0], disp[1], disp[2]);
			}
			return vec3d(0.0, 0.0, 0.0);
		});
		return true;
	}

private:
	const DisplacementContainer& m_data;
};

// -----------------------------------------------------------------------------'

class VFMPlotDeformationGradient : public FEPlotDomainData
{
public:
	VFMPlotDeformationGradient(FEModel* fem, const DeformationGradientField& field)
		: FEPlotDomainData(fem, PLT_MAT3F, FMT_ITEM), m_field(field)
	{
		SetUnits(UNIT_NONE);
	}

	bool Save(FEDomain& dom, FEDataStream& a) override
	{
		if (dom.Class() != FE_DOMAIN_SOLID)
		{
			// Fill with identity matrices for non-solid domains
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
			const GaussPointDeformation* gpF = m_field.Find(el.GetID());

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
	const DeformationGradientField& m_field;
};

} // namespace

// -----------------------------------------------------------------------------'

bool ExportVFMKinematics(const std::string& filePath,
	FEModel& fem,
	const DisplacementContainer& measured,
	const DisplacementContainer& virt,
	const DeformationGradientField& defGrad,
	std::string& error)
{
	VFMPlotFile plt(&fem);

	auto* measuredField = new VFMPlotMeasuredDisplacement(&fem, measured);
	if (!plt.AddVariable(measuredField, "vfm_measured_disp"))
	{
		delete measuredField;
		error = "Failed to register measured displacement field.";
		return false;
	}

	auto* virtualField = new VFMPlotVirtualDisplacement(&fem, virt);
	if (!plt.AddVariable(virtualField, "vfm_virtual_disp"))
	{
		delete virtualField;
		error = "Failed to register virtual displacement field.";
		return false;
	}

	auto* defGradField = new VFMPlotDeformationGradient(&fem, defGrad);
	if (!plt.AddVariable(defGradField, "vfm_deformation_gradient"))
	{
		delete defGradField;
		error = "Failed to register deformation gradient field.";
		return false;
	}

	plt.SetSoftwareString("FEBio VFM plug-in");

	if (!plt.Open(filePath.c_str()))
	{
		error = "Unable to create plot file: " + filePath;
		return false;
	}

	if (!plt.Write(0.0f))
	{
		error = "Failed to write plot state.";
		plt.Close();
		return false;
	}

	plt.Close();
	return true;
}
