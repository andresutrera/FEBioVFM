#include "FE/shape_provider_febio.hpp"
#include <FECore/FESolidDomain.h>
#include <FECore/FESolidElement.h>

void FEBioShapeProvider::gradN(size_t e, size_t g, std::vector<vec3d>& dNdx0) const {
  FESolidDomain* dom = _c.elemRef[e].dom;
  const int k = _c.elemRef[e].local;
  FESolidElement& el = static_cast<FESolidElement&>(dom->ElementRef(k));

  const int nen = el.Nodes();
  dNdx0.resize(nen);

  double Ji[3][3]; dom->invjac0(el, Ji, (int)g);
  double* Gr = el.Gr((int)g);
  double* Gs = el.Gs((int)g);
  double* Gt = el.Gt((int)g);

  for (int a=0; a<nen; ++a) {
    const double GX = Ji[0][0]*Gr[a] + Ji[1][0]*Gs[a] + Ji[2][0]*Gt[a];
    const double GY = Ji[0][1]*Gr[a] + Ji[1][1]*Gs[a] + Ji[2][1]*Gt[a];
    const double GZ = Ji[0][2]*Gr[a] + Ji[1][2]*Gs[a] + Ji[2][2]*Gt[a];
    dNdx0[a] = vec3d(GX,GY,GZ);
  }
}
