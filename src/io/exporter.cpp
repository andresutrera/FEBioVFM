#include "io/exporter.hpp"
#include "optimization/vfm_problem.hpp"
#include "domain/vfm_displacements.hpp"
#include "domain/vfm_tensors.hpp"
#include <FEBioPlot/FEBioPlotFile.h>
#include <FECore/writeplot.h>
#include <FECore/FEMesh.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FESolidElement.h>
#include <FECore/units.h>
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <algorithm>
#include <utility>

namespace
{

  class VFMPlotFile : public FEBioPlotFile
  {
  public:
    explicit VFMPlotFile(FEModel *fem) : FEBioPlotFile(fem) {}
    using PlotFile::AddVariable;
  };

  class ProblemLogger
  {
  public:
    explicit ProblemLogger(FEModel &model) : fem(&model) {}
    FEModel *GetFEModel() const { return fem; }

    template <typename... Args>
    void info(const char *fmt, Args &&...args) const
    {
      feLog(fmt, std::forward<Args>(args)...);
    }

  private:
    FEModel *fem = nullptr;
  };

  class DisplacementPlot : public FEPlotNodeData
  {
  public:
    DisplacementPlot(FEModel *fem, const MeshDims &dims)
        : FEPlotNodeData(fem, PLT_VEC3F, FMT_NODE), m_dims(&dims)
    {
      SetUnits(UNIT_LENGTH);
    }

    void setField(const NodalField<vec3d> *field) { m_field = field; }

    bool Save(FEMesh &mesh, FEDataStream &a) override
    {
      writeNodalValues<vec3d>(mesh, a, [&](const FENode &node) -> vec3d
                              {
      if (!m_field)
        return vec3d(0.0, 0.0, 0.0);
      auto it = m_dims->nodeId2idx.find(node.GetID());
      if (it == m_dims->nodeId2idx.end())
        return vec3d(0.0, 0.0, 0.0);
      std::size_t idx = it->second;
      if (idx >= m_field->size())
        return vec3d(0.0, 0.0, 0.0);
      return m_field->getNode(idx); });
      return true;
    }

  private:
    const MeshDims *m_dims = nullptr;
    const NodalField<vec3d> *m_field = nullptr;
  };

  class DeformationGradientPlot : public FEPlotDomainData
  {
  public:
    DeformationGradientPlot(FEModel *fem, const MeshDims &dims)
        : FEPlotDomainData(fem, PLT_MAT3F, FMT_ITEM), m_dims(&dims)
    {
      SetUnits(UNIT_NONE);
    }

    void setField(const RaggedElemField<mat3d> *field) { m_field = field; }

    bool Save(FEDomain &dom, FEDataStream &a) override
    {
      if (dom.Class() != FE_DOMAIN_SOLID)
      {
        for (int i = 0; i < dom.Elements(); ++i)
          a << mat3d::identity();
        return true;
      }

      FESolidDomain &sd = static_cast<FESolidDomain &>(dom);
      for (int i = 0; i < sd.Elements(); ++i)
      {
        FESolidElement &el = sd.Element(i);
        mat3d avg = mat3d::identity();
        if (m_field)
        {
          auto it = m_dims->elemId2idx.find(el.GetID());
          if (it != m_dims->elemId2idx.end())
          {
            std::size_t eidx = it->second;
            const std::size_t nint = m_field->nGauss((ElemIdx)eidx);
            if (nint > 0)
            {
              mat3d sum;
              sum.zero();
              for (std::size_t g = 0; g < nint; ++g)
                sum += m_field->getElemGP((ElemIdx)eidx, (GPIdx)g);
              sum /= static_cast<double>(nint);
              avg = sum;
            }
          }
        }
        a << avg;
      }
      return true;
    }

  private:
    const MeshDims *m_dims = nullptr;
    const RaggedElemField<mat3d> *m_field = nullptr;
  };

  class CauchyStressPlot : public FEPlotDomainData
  {
  public:
    CauchyStressPlot(FEModel *fem, const MeshDims &dims)
        : FEPlotDomainData(fem, PLT_MAT3FS, FMT_ITEM), m_dims(&dims)
    {
      SetUnits(UNIT_PRESSURE);
    }

    void setField(const RaggedElemField<mat3d> *sigma) { m_sigma = sigma; }

    bool Save(FEDomain &dom, FEDataStream &a) override
    {
      mat3ds zero;
      zero.zero();
      if (dom.Class() != FE_DOMAIN_SOLID)
      {
        for (int i = 0; i < dom.Elements(); ++i)
          a << zero;
        return true;
      }

      FESolidDomain &sd = static_cast<FESolidDomain &>(dom);
      for (int i = 0; i < sd.Elements(); ++i)
      {
        mat3ds avg;
        avg.zero();
        if (m_sigma)
        {
          auto it = m_dims->elemId2idx.find(sd.Element(i).GetID());
          if (it != m_dims->elemId2idx.end())
          {
            std::size_t eidx = it->second;
            const std::size_t nint = m_sigma->nGauss((ElemIdx)eidx);
            if (nint > 0)
            {
              for (std::size_t g = 0; g < nint; ++g)
              {
                const mat3d &m = m_sigma->getElemGP((ElemIdx)eidx, (GPIdx)g);
                mat3ds s = m.sym();
                avg += s;
              }
              avg /= static_cast<double>(nint);
            }
          }
        }
        a << avg;
      }
      return true;
    }

  private:
    const MeshDims *m_dims = nullptr;
    const RaggedElemField<mat3d> *m_sigma = nullptr;
  };

  class FirstPiolaPlot : public FEPlotDomainData
  {
  public:
    FirstPiolaPlot(FEModel *fem, const MeshDims &dims)
        : FEPlotDomainData(fem, PLT_MAT3F, FMT_ITEM), m_dims(&dims)
    {
      SetUnits(UNIT_PRESSURE);
    }

    void setField(const RaggedElemField<mat3d> *piola) { m_piola = piola; }

    bool Save(FEDomain &dom, FEDataStream &a) override
    {
      mat3d zero;
      zero.zero();
      if (dom.Class() != FE_DOMAIN_SOLID)
      {
        for (int i = 0; i < dom.Elements(); ++i)
          a << zero;
        return true;
      }

      FESolidDomain &sd = static_cast<FESolidDomain &>(dom);
      for (int i = 0; i < sd.Elements(); ++i)
      {
        mat3d avg;
        avg.zero();
        if (m_piola)
        {
          auto it = m_dims->elemId2idx.find(sd.Element(i).GetID());
          if (it != m_dims->elemId2idx.end())
          {
            std::size_t eidx = it->second;
            const std::size_t nint = m_piola->nGauss((ElemIdx)eidx);
            if (nint > 0)
            {
              for (std::size_t g = 0; g < nint; ++g)
                avg += m_piola->getElemGP((ElemIdx)eidx, (GPIdx)g);
              avg /= static_cast<double>(nint);
            }
          }
        }
        a << avg;
      }
      return true;
    }

  private:
    const MeshDims *m_dims = nullptr;
    const RaggedElemField<mat3d> *m_piola = nullptr;
  };

  struct VirtualFieldPlots
  {
    DisplacementPlot *disp = nullptr;
    DeformationGradientPlot *def = nullptr;
    bool dispSingle = false;
    bool defSingle = false;
  };

  std::string vf_name(const std::string &base, const std::string &id, std::size_t idx, std::size_t total)
  {
    if (!id.empty())
      return base + " " + id;
    if (total > 1)
      return base + " #" + std::to_string(idx);
    return base;
  }

} // namespace

bool export_vfm_results(const VFMProblem &problem,
                        const std::string &filePath,
                        std::string &err)
{
  err.clear();
  if (problem.fem == nullptr)
  {
    err = "VFM problem not initialized.";
    return false;
  }

  const std::size_t measTimes = problem.state.measured.series.nTimes();
  const std::size_t defTimes = problem.state.def.nTimes();
  const std::size_t stressTimes = problem.state.stresses.nTimes();

  std::size_t maxTimes = std::max(measTimes, std::max(defTimes, stressTimes));
  for (std::size_t v = 0; v < problem.state.virtuals.nVF(); ++v)
    maxTimes = std::max<std::size_t>(maxTimes, problem.state.virtuals.getVF((VFIdx)v).nTimes());
  for (std::size_t v = 0; v < problem.state.vdef.nVF(); ++v)
    maxTimes = std::max<std::size_t>(maxTimes, problem.state.vdef.getVF((VFIdx)v).nTimes());

  if (maxTimes == 0)
  {
    err = "No frames available for export.";
    return false;
  }

  ProblemLogger logger(*problem.fem);

  VFMPlotFile plot(problem.fem);
  plot.SetSoftwareString("FEBio VFM exporter");

  // Measured displacement
  auto *measDisp = new DisplacementPlot(problem.fem, problem.dims);
  if (!plot.AddVariable(measDisp, "displacement"))
  {
    err = "Failed to register measured displacement.";
    delete measDisp;
    return false;
  }

  // Measured deformation gradient
  auto *measDef = new DeformationGradientPlot(problem.fem, problem.dims);
  if (!plot.AddVariable(measDef, "measured deformation gradient"))
  {
    err = "Failed to register measured deformation gradient.";
    delete measDef;
    return false;
  }

  // Stress plots
  auto *cauchyPlot = new CauchyStressPlot(problem.fem, problem.dims);
  if (!plot.AddVariable(cauchyPlot, "cauchy stress"))
  {
    err = "Failed to register Cauchy stress.";
    delete cauchyPlot;
    return false;
  }

  auto *piolaPlot = new FirstPiolaPlot(problem.fem, problem.dims);
  if (!plot.AddVariable(piolaPlot, "first piola stress"))
  {
    err = "Failed to register First Piola stress.";
    delete piolaPlot;
    return false;
  }

  // Virtual field plots
  std::vector<VirtualFieldPlots> vfPlots;
  const std::size_t nVF = problem.state.virtuals.nVF();
  vfPlots.reserve(nVF);
  for (std::size_t v = 0; v < nVF; ++v)
  {
    const auto &vf = problem.state.virtuals.getVF((VFIdx)v);
    const auto &vfDefSeries = problem.state.vdef.getVF((VFIdx)v);
    if (vf.nTimes() == 0)
    {
      err = "virtual displacement field has no time steps.";
      return false;
    }
    if (vfDefSeries.nTimes() == 0)
    {
      err = "virtual deformation field has no time steps.";
      return false;
    }
    std::string dispName = vf_name("virtual displacement", "", v, nVF);
    auto *vd = new DisplacementPlot(problem.fem, problem.dims);
    if (!plot.AddVariable(vd, dispName.c_str()))
    {
      err = "Failed to register virtual displacement variable.";
      delete vd;
      return false;
    }

    std::string defName = vf_name("virtual deformation gradient", "", v, nVF);
    auto *vg = new DeformationGradientPlot(problem.fem, problem.dims);
    if (!plot.AddVariable(vg, defName.c_str()))
    {
      err = "Failed to register virtual deformation gradient variable.";
      delete vg;
      return false;
    }

    const bool dispSingle = (vf.nTimes() == 1);
    const bool defSingle = (vfDefSeries.nTimes() == 1);
    vfPlots.push_back({vd, vg, dispSingle, defSingle});
  }

  if (!plot.Open(filePath.c_str()))
  {
    err = "Unable to create plot file: " + filePath;
    return false;
  }

  for (std::size_t t = 0; t < maxTimes; ++t)
  {
    // measured displacement
    const NodalField<vec3d> *dispField = nullptr;
    if (t < measTimes)
      dispField = &problem.state.measured.series.getTime((TimeIdx)t).u;
    measDisp->setField(dispField);

    // measured deformation gradient
    const RaggedElemField<mat3d> *defField = nullptr;
    if (t < defTimes)
      defField = &problem.state.def.series.getTime((TimeIdx)t).F;
    measDef->setField(defField);

    // stresses
    const RaggedElemField<mat3d> *sigmaField = nullptr;
    const RaggedElemField<mat3d> *piolaFieldPtr = nullptr;
    if (t < stressTimes)
    {
      const auto &frame = problem.state.stresses.series.getTime((TimeIdx)t);
      sigmaField = &frame.sigma;
      piolaFieldPtr = &frame.P;
    }
    cauchyPlot->setField(sigmaField);
    piolaPlot->setField(piolaFieldPtr);

    // virtual fields
    for (std::size_t v = 0; v < vfPlots.size(); ++v)
    {
      const auto &vfDispSeries = problem.state.virtuals.getVF((VFIdx)v);
      const auto &vfDefSeries = problem.state.vdef.getVF((VFIdx)v);

      const NodalField<vec3d> *vdisp = nullptr;
      if (vfPlots[v].dispSingle)
      {
        if (vfDispSeries.nTimes() > 0)
          vdisp = &vfDispSeries.getTime((TimeIdx)0).u;
      }
      else if (t < static_cast<std::size_t>(vfDispSeries.nTimes()))
      {
        vdisp = &vfDispSeries.getTime((TimeIdx)t).u;
      }
      vfPlots[v].disp->setField(vdisp);

      const RaggedElemField<mat3d> *vdef = nullptr;
      if (vfPlots[v].defSingle)
      {
        if (vfDefSeries.nTimes() > 0)
          vdef = &vfDefSeries.getTime((TimeIdx)0).F;
      }
      else if (t < static_cast<std::size_t>(vfDefSeries.nTimes()))
      {
        vdef = &vfDefSeries.getTime((TimeIdx)t).F;
      }
      vfPlots[v].def->setField(vdef);
    }

    const double timeValue = static_cast<double>(t);
    if (!plot.Write(static_cast<float>(timeValue)))
    {
      err = "Failed to write plot state.";
      plot.Close();
      return false;
    }
  }

  plot.Close();
  logger.info("Exported XPLT results to %s\n", filePath.c_str());
  return true;
}
