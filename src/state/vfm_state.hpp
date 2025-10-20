// src/state/vfm_state.hpp
#pragma once
#include <string>
#include <vector>

// displacements
#include "domain/vfm_displacements.hpp"   // MeasuredData, VirtualFields, MeasuredLoad
// tensors
#include "domain/vfm_tensors.hpp"

struct VFMParamSpec {
  std::string name;
  double init  = 0.0;
  double lo    = 0.0;
  double hi    = 0.0;
  double scale = 1.0;
};

struct VFMParam {
  VFMParamSpec spec;   // XML spec: name, init, bounds, scale
  double value = 0.0;  // current physical value used in the model
};

struct VFMState {
    // ---------------- inputs ----------------
    MeasuredData   measured;    // u(t,i)
    VirtualFields  virtuals;    // u_v(v,t,i)
    MeasuredLoad   loads;       // F(t,i)  (optional, can be empty)

    // ------------- derived tensors ----------
    Deformations         def;    // F(t,e,g)
    VirtualDeformations  vdef;   // F_v(v,t,e,g)
    Stresses             stresses; // σ,P(t,e,g) (optional for later)

    std::vector<VFMParam> params;


    // utilities
    void clear() { *this = VFMState{}; }

    // size tensors once you know quadrature shape and VF count
    void configure_tensors(const std::vector<std::size_t>& gpPerElem, std::size_t nVF) {
        def.setElemShape(gpPerElem);
        vdef.setElemShape(gpPerElem);
        vdef.resizeVF(nVF);
        stresses.setElemShape(gpPerElem);
    }

    // create time frames to mirror current displacement timelines
    // call AFTER configure_tensors
    void mirror_frames_from_displacements() {
        // measured → def
        for (std::size_t k=0; k<measured.series.nTimes(); ++k) (void)def.addTime();
        // virtuals → vdef
        for (VFIdx v=0; v<(VFIdx)virtuals.nVF(); ++v)
            for (std::size_t k=0; k<virtuals.getVF(v).nTimes(); ++k)
                (void)vdef.addTime(v);
    }

    std::vector<double> pack_theta() const {
        std::vector<double> t; t.reserve(params.size());
        for (const auto& q : params) t.push_back((q.value - q.spec.init) / q.spec.scale);
        return t;
    }
    void unpack_theta(const std::vector<double>& t) {
        for (size_t i=0;i<params.size();++i)
        params[i].value = params[i].spec.init + params[i].spec.scale * t[i];
    }
};
