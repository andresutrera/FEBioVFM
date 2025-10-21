// fe/material_provider_febio.cpp
#include "FE/material_provider_febio.hpp"
#include <FECore/FESolidDomain.h>
#include <FECore/FESolidElement.h>
#include <FEBioMech/FESolidMaterial.h>
#include <FEBioMech/FEUncoupledMaterial.h>
#include <FECore/FECoreKernel.h>
#include <FEBioMech/FEElasticMaterialPoint.h>
#include <iostream>
bool FEBioMaterialProvider::evalCauchy(size_t e, size_t g, const mat3d &F, mat3ds &sigma, std::string &err) const
{
  auto *dom = _c.elemRef[e].dom;
  const int k = _c.elemRef[e].local;
  FESolidElement &el = static_cast<FESolidElement &>(dom->ElementRef(k));

  auto *base = dynamic_cast<FESolidMaterial *>(dom->GetMaterial());
  if (!base)
  {
    err = "domain has no FESolidMaterial";
    return false;
  }

  FEMaterialPoint *mp0 = el.GetMaterialPoint((int)g);
  if (!mp0)
  {
    err = "null material point";
    return false;
  }

  std::unique_ptr<FEMaterialPoint> mp(mp0->Copy());
  if (!mp)
  {
    err = "failed to clone material point";
    return false;
  }

  mp->m_elem = &el;
  mp->m_index = (int)g;

  auto *ep = mp->ExtractData<FEElasticMaterialPoint>();
  if (!ep)
  {
    err = "no FEElasticMaterialPoint";
    return false;
  }

  ep->m_F = F;
  ep->m_J = F.det();
  ep->m_s.zero();
  ep->m_L.zero();
  ep->m_v = ep->m_a = vec3d(0, 0, 0);
  ep->m_Wt = ep->m_Wp = 0.0;

  if (auto *unc = dynamic_cast<FEUncoupledMaterial *>(base))
  {
    // Mirror legacy VFMStress path where deviatoric stress is adjusted by a
    // pressure-like term to recover the full Cauchy stress.
    // This assumes sigma zz = 0!!!!!!
    // TODO: Fix this to dont loose generality.
    mat3ds dev = unc->DevStress(*mp);
    const double p = dev.zz();
    sigma = dev - mat3dd(p);
  }
  else
  {
    sigma = base->Stress(*mp);
  }

  return true;
}
