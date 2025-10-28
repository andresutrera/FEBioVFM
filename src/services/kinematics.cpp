#include "services/kinematics.hpp"
#include <cmath>

static inline mat3d F_at(size_t e, size_t g, const IShapeProvider &sp, const NodalField<vec3d> &ut)
{
    const auto &en = sp.elemNodes(e);
    std::vector<vec3d> dN;
    sp.gradN(e, g, dN);
    mat3d F;
    F.unit();
    for (size_t a = 0; a < en.size(); ++a)
    {
        const vec3d u = ut.getNode(en[a]), G = dN[a];
        F[0][0] += u.x * G.x;
        F[0][1] += u.x * G.y;
        F[0][2] += u.x * G.z;
        F[1][0] += u.y * G.x;
        F[1][1] += u.y * G.y;
        F[1][2] += u.y * G.z;
        F[2][0] += u.z * G.x;
        F[2][1] += u.z * G.y;
        F[2][2] += u.z * G.z;
    }
    return F;
}
static inline bool set_F_checked(Deformations &Fo, TimeIdx t, size_t e, size_t g,
                                 const mat3d &F, bool checkDet, std::string &err)
{
    if (checkDet && F.det() <= 0.0)
    {
        err = "non-positive det(F)";
        return false;
    }
    Fo.setF(t, e, g, F);
    return true;
}

static inline bool set_Fv_checked(VirtualDeformations &Fv, VFIdx v, TimeIdx t, size_t e, size_t g,
                                  const mat3d &F, bool checkDet, std::string &err)
{
    if (checkDet && F.det() <= 0.0)
    {
        err = "non-positive det(F)";
        return false;
    }
    Fv.setF(v, t, e, g, F);
    return true;
}

bool Kinematics::compute_measured(const MeshQuad &q, const IShapeProvider &shp,
                                  const MeasuredData &U, Deformations &Fo,
                                  bool planeDeformation,
                                  bool checkDet, std::string &err)
{
    for (TimeIdx t = 0; t < (TimeIdx)U.series.nTimes(); ++t)
    {
        const auto &ut = U.series.getTime(t).u;
        for (size_t e = 0; e < q.gpPerElem.size(); ++e)
        {
            for (size_t g = 0; g < q.gpPerElem[e]; ++g)
            {
                mat3d F = F_at(e, g, shp, ut);
                if (planeDeformation)
                {
                    F[0][2] = 0.0;
                    F[1][2] = 0.0;
                    F[2][0] = 0.0;
                    F[2][1] = 0.0;
                    const double denom = F[0][0] * F[1][1];
                    F[2][2] = 1.0 / denom;
                }
                if (!set_F_checked(Fo, t, e, g, F, checkDet, err))
                    return false;
            }
        }
    }
    return true;
}

bool Kinematics::compute_virtuals(const MeshQuad &q, const IShapeProvider &shp,
                                  const VirtualFields &UV, VirtualDeformations &FV,
                                  bool checkDet, std::string &err)
{
    for (VFIdx v = 0; v < (VFIdx)UV.nVF(); ++v)
    {
        const auto &ts = UV.getVF(v);
        for (TimeIdx t = 0; t < (TimeIdx)ts.nTimes(); ++t)
        {
            const auto &ut = ts.getTime(t).u;
            for (size_t e = 0; e < q.gpPerElem.size(); ++e)
            {
                for (size_t g = 0; g < q.gpPerElem[e]; ++g)
                {
                    const mat3d F = F_at(e, g, shp, ut);
                    if (!set_Fv_checked(FV, v, t, e, g, F, checkDet, err))
                        return false;
                }
            }
        }
    }
    return true;
}
