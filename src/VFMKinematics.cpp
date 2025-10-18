#include "VFMKinematics.h"

#include <array>
#include <vector>

#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FESolidElement.h>
#include <FECore/log.h>
#include <FECore/mat3d.h>

#include "DisplacementContainer.h"
#include "DeformationGradientField.h"

namespace {

/*
 * The deformation gradient is reconstructed directly from the nodal displacement
 * field without mutating the FEBio mesh.  Given displacements \f$\mathbf{u}_i\f$
 * defined at the element nodes and the reference gradients of the shape
 * functions \f$\nabla_X N_i\f$, the classical Total Lagrangian relation
 *
 * \f[
 *   \mathbf{F} = \mathbf{I} + \sum_{i=1}^{n_{\text{eln}}} \mathbf{u}_i \otimes
 *   \nabla_X N_i
 * \f]
 *
 * is applied at every Gauss point.  The reference gradients are obtained from
 * the inverse Jacobian \f$\mathbf{J}_0^{-1}\f$ (evaluated in the reference
 * configuration via \c invjac0).  When all displacements are zero the summation
 * vanishes and \f$\mathbf{F} = \mathbf{I}\f$, which serves as a convenient
 * sanity check.
 */
bool ComputeDefGrad(
	FESolidDomain& domain,
	FEModel& fem,
	FESolidElement& el,
	const std::vector<vec3d>& u,
	int n,
	mat3d& F,
	std::string& errorMessage)
{
	double Ji[3][3];
	domain.invjac0(el, Ji, n);

	double* Grn = el.Gr(n);
	double* Gsn = el.Gs(n);
	double* Gtn = el.Gt(n);

	F.zero();

	const int neln = el.Nodes();
	for (int i = 0; i < neln; ++i)
	{
		const double Gri = Grn[i];
		const double Gsi = Gsn[i];
		const double Gti = Gtn[i];

		const double GX = Ji[0][0] * Gri + Ji[1][0] * Gsi + Ji[2][0] * Gti;
		const double GY = Ji[0][1] * Gri + Ji[1][1] * Gsi + Ji[2][1] * Gti;
		const double GZ = Ji[0][2] * Gri + Ji[1][2] * Gsi + Ji[2][2] * Gti;

		const vec3d& ui = u[i];

		F[0][0] += GX * ui.x; F[0][1] += GY * ui.x; F[0][2] += GZ * ui.x;
		F[1][0] += GX * ui.y; F[1][1] += GY * ui.y; F[1][2] += GZ * ui.y;
		F[2][0] += GX * ui.z; F[2][1] += GY * ui.z; F[2][2] += GZ * ui.z;
	}

	F[0][0] += 1.0;
	F[1][1] += 1.0;
	F[2][2] += 1.0;

	const double detF = F.det();

	feLogDebugEx(&fem,
		"F(elem %d, gp %d) = [[%.5f %.5f %.5f] [%.5f %.5f %.5f] [%.5f %.5f %.5f]], det(F) = %.5f",
		el.GetID(), n,
		F[0][0], F[0][1], F[0][2],
		F[1][0], F[1][1], F[1][2],
		F[2][0], F[2][1], F[2][2],
		detF);

	if (detF <= 0.0)
	{
		errorMessage = "Computed deformation gradient has non-positive determinant.";
		return false;
	}

	return true;
}

} // namespace

bool VFMKinematics::ComputeDeformationGradients(FEModel& fem,
	const DisplacementContainer& displacements,
	DeformationGradientField& outField,
	std::string& errorMessage)
{
	outField.Clear();

	FEMesh& mesh = fem.GetMesh();

	// Look all domains
	const int domainCount = mesh.Domains();
	for (int domIdx = 0; domIdx < domainCount; ++domIdx)
	{

		FEDomain& domain = mesh.Domain(domIdx);
		auto* solidDomain = dynamic_cast<FESolidDomain*>(&domain);
		if (solidDomain == nullptr) continue;

		// Look domain elements
		const int elemCount = solidDomain->Elements();
		for (int elemIdx = 0; elemIdx < elemCount; ++elemIdx)
		{
			FESolidElement& el = static_cast<FESolidElement&>(solidDomain->ElementRef(elemIdx));
			const int neln = el.Nodes();

			std::vector<vec3d> displacementsVec(neln);
			// Look element nodes
			for (int i = 0; i < neln; ++i)
			{
				const int nodeIndex = el.m_node[i];
				const int nodeId = mesh.Node(nodeIndex).GetID();
				std::array<double, 3> disp{};
				if (!displacements.TryGet(nodeId, disp))
				{
					errorMessage = "Missing displacement entry for node " + std::to_string(nodeId) + ".";
					return false;
				}

				const vec3d u(disp[0], disp[1], disp[2]);
				displacementsVec[i] = u;
			}

			GaussPointDeformation gpData;
			gpData.elementId = el.GetID();
			const int nint = el.GaussPoints();
			gpData.gradients.resize(nint);


			feLogDebugEx(&fem,"Computed deformation gradient field:");
			// For every gauss point
			for (int n = 0; n < nint; ++n)
			{
				mat3d F;
				if (!ComputeDefGrad(*solidDomain, fem, el, displacementsVec, n, F, errorMessage))
				{
					errorMessage += " Element ID: " + std::to_string(gpData.elementId) + ", integration point: " + std::to_string(n);
					return false;
				}
				gpData.gradients[n] = F;
			}

			outField.Add(std::move(gpData));
		}
	}

	return true;
}
