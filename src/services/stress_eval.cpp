#include "services/stress_eval.hpp"
bool StressEval::cauchy(const Deformations& F, Stresses& out, const IMaterialProvider& mat, std::string& err){
  if (out.nTimes()!=F.nTimes()) { // assume frames already created to match
    for (size_t k=out.nTimes(); k<F.nTimes(); ++k) (void)out.addTime();
  }
  for (TimeIdx t=0; t<(TimeIdx)F.nTimes(); ++t){
    const size_t ne = F.nElements(t);
    for (size_t e=0; e<ne; ++e){
      const size_t ng = F.nGauss(t,e);
      for (size_t g=0; g<ng; ++g){
        const mat3d& Fin = F.crefF(t,e,g);
        mat3ds s;
        if (!mat.evalCauchy(e,g, Fin, s, err)) return false;
        out.setSigma(t,e,g, mat3d( // store as mat3d or keep mat3ds via an overload
          s.xx(), s.xy(), s.xz(),
          s.xy(), s.yy(), s.yz(),
          s.xz(), s.yz(), s.zz()));
      }
    }
  }
  return true;
}

bool StressEval::first_piola(const Deformations& F, const Stresses& S, Stresses& outP, std::string& err){
  if (outP.nTimes()!=F.nTimes()) {
    for (size_t k=outP.nTimes(); k<F.nTimes(); ++k) (void)outP.addTime();
  }
  for (TimeIdx t=0; t<(TimeIdx)F.nTimes(); ++t){
    const size_t ne = F.nElements(t);
    for (size_t e=0; e<ne; ++e){
      const size_t ng = F.nGauss(t,e);
      for (size_t g=0; g<ng; ++g){
        const mat3d& Fin = F.crefF(t,e,g);
        const double J = Fin.det();
        if (J<=0) { err = "non-positive detF"; return false; }
        const mat3d FinvT = Fin.inverse().transpose();

        const mat3d Sm = S.crefSigma(t,e,g); // stored as full
        mat3d P = Sm * FinvT; P *= J;
        outP.setP(t,e,g, P);
      }
    }
  }
  return true;
}
