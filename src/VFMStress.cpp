#include "VFMStress.h"

#include "DeformationGradientField.h"
#include "StressField.h"

#include <memory>
#include <string>

#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FESolidElement.h>
#include <FECore/mat3d.h>
#include <FECore/log.h>

#include <FEBioMech/FESolidMaterial.h>
#include <FEBioMech/FEUncoupledMaterial.h>
#include <FEBioMech/FEElasticMaterialPoint.h>

bool VFMStress::ComputeCauchyStress(FEModel& fem,
	const DeformationGradientField& defField,
	StressField& outField,
	std::string& errorMessage)
{
	outField.Clear();

	FEMesh& mesh = fem.GetMesh();

	const int domainCount = mesh.Domains();
	for (int domIdx = 0; domIdx < domainCount; ++domIdx)
	{
		FEDomain& domain = mesh.Domain(domIdx);
		auto* solidDomain = dynamic_cast<FESolidDomain*>(&domain);
		if (solidDomain == nullptr) continue;

		auto* solidMaterial = dynamic_cast<FESolidMaterial*>(solidDomain->GetMaterial());
		if (solidMaterial == nullptr)
		{
			errorMessage = "Encountered a solid domain without a compatible solid material instance.";
			return false;
		}

		const int elemCount = solidDomain->Elements();
		for (int elemIdx = 0; elemIdx < elemCount; ++elemIdx)
		{
			FESolidElement& el = static_cast<FESolidElement&>(solidDomain->ElementRef(elemIdx));
			const int elemId = el.GetID();

			const GaussPointDeformation* gpDef = defField.Find(elemId);
			if (gpDef == nullptr)
			{
				errorMessage = "Missing deformation gradient data for element " + std::to_string(elemId) + ".";
				return false;
			}

			const int gaussPointCount = el.GaussPoints();
			if ((int)gpDef->gradients.size() != gaussPointCount)
			{
				errorMessage = "Mismatch between stored deformation gradients and element integration points for element " + std::to_string(elemId) + ".";
				return false;
			}

			GaussPointStress gpStress;
			gpStress.elementId = elemId;
			gpStress.stresses.resize(gaussPointCount);

			for (int n = 0; n < gaussPointCount; ++n)
			{
				const mat3d& F = gpDef->gradients[n];

				FEMaterialPoint* originalPoint = el.GetMaterialPoint(n);
				if (originalPoint == nullptr)
				{
					errorMessage = "Element " + std::to_string(elemId) + " lacks material point data at integration point " + std::to_string(n) + ".";
					return false;
				}

				std::unique_ptr<FEMaterialPoint> mpClone(originalPoint->Copy());
				if (!mpClone)
				{
					errorMessage = "Failed to clone material point state for element " + std::to_string(elemId) + ".";
					return false;
				}

				mpClone->m_elem = &el;
				mpClone->m_index = n;

				FEElasticMaterialPoint* elasticPoint = mpClone->ExtractData<FEElasticMaterialPoint>();
				if (elasticPoint == nullptr)
				{
					errorMessage = "Material in element " + std::to_string(elemId) + " does not expose elastic material point data.";
					return false;
				}

				elasticPoint->m_F = F;
				elasticPoint->m_J = F.det();
				elasticPoint->m_s.zero();
				elasticPoint->m_v = elasticPoint->m_a = vec3d(0.0, 0.0, 0.0);
				elasticPoint->m_L.zero();
				elasticPoint->m_Wt = elasticPoint->m_Wp = 0.0;

				mat3ds sigma;
				if (auto* uncoupled = dynamic_cast<FEUncoupledMaterial*>(solidMaterial))
				{
					mat3ds dev = uncoupled->DevStress(*mpClone);
					const double p = dev.zz();
					sigma = dev - mat3dd(p);
				}
				else
				{
					sigma = solidMaterial->Stress(*mpClone);
				}
				gpStress.stresses[n] = sigma;
			}

			outField.Add(std::move(gpStress));
		}
	}

	return true;
}
