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

/**
 * \brief Reconstruct the deformation gradient purely from nodal displacements.
 *
 * The computation follows the Total Lagrangian relation
 * \f[
 *   \mathbf{F} = \mathbf{I} + \sum_{i=1}^{n_{\mathrm{eln}}} \mathbf{u}_i \otimes \nabla_X N_i,
 * \f]
 * where \f$\mathbf{u}_i\f$ denotes the nodal displacement vectors and
 * \f$\nabla_X N_i\f$ the gradients of the shape functions with respect to the
 * reference configuration. The gradients are assembled using the inverse
 * reference Jacobian \f$\mathbf{J}_0^{-1}\f$ obtained via
 * `FESolidDomain::invjac0`. When all displacements are zero the summation
 * vanishes and \f$\mathbf{F} = \mathbf{I}\f`, providing an immediate
 * consistency check.
 *
 * @param domain Solid domain that supplies geometric data and reference Jacobians.
 * @param fem Owning FEBio model used for debug logging.
 * @param el Element whose deformation gradient is reconstructed.
 * @param u Nodal displacement vectors expressed in the global frame.
 * @param n Index of the Gauss point inside \p el.
 * @param[out] F Resulting deformation gradient.
 * @param[out] errorMessage Filled with diagnostic text when the determinant is non-positive.
 * @return `true` on success, `false` when the determinant is non-positive or a debug failure occurs.
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
            const int neln = el.Nodes();

            std::vector<vec3d> displacementsVec(neln);
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

                displacementsVec[i] = vec3d(disp[0], disp[1], disp[2]);
            }

            GaussPointDeformation gpData;
            gpData.elementId = el.GetID();
            const int nint = el.GaussPoints();
            gpData.gradients.resize(nint);

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
