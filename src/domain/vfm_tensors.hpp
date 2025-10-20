// vfm_tensor_series.hpp
#pragma once
#include <vector>
#include <cstddef>
#include <cassert>
#include <FECore/mat3d.h> // kept local to the tensors side
#include "domain/vfm_core_series.hpp"
// ------------------------------ indices -------------------------------------
// Shared small POD aliases. Keep these consistent across modules.
using ElemIdx = std::size_t; // element index in [0..nElem)
using GPIdx   = std::size_t; // gauss-point index local to element e
using TimeIdx = int;         // time/frame index in [0..nTimes)
using VFIdx   = int;         // virtual-field index in [0..nVF)

// ========================= RAGGED ELEM×GP STORAGE ============================
// RaggedElemField<T>
// Purpose:
//   Compact CSR-like storage for element-wise arrays with NON-UNIFORM
//   gauss-point counts. Supports:
//     • one-shot shape: setElemShape(vector<nGP_e>)
//     • two-pass shape: prepare(nElem) → setGaussCount(e,ng) → finalize()
// Guarantees:
//   After finalize/setElemShape, storage is fixed-size. get/set are O(1).
// Notes:
//   T is the stored tensor type (e.g., mat3d).
// ============================================================================
template <typename T>
class RaggedElemField {
public:
    // -------- shape build --------
    void setElemShape(const std::vector<std::size_t>& nGPperElem) { _build(nGPperElem); }

    void prepare(std::size_t nElem) {
        _nGP.assign(nElem, 0);
        _ofs.assign(nElem + 1, 0);
        _data.clear();
        _final = false;
    }
    void setGaussCount(ElemIdx e, std::size_t ng) {
        assert(!_final && e < _nGP.size());
        _nGP[e] = ng;
    }
    void finalize() {
        assert(!_final);
        for (std::size_t e = 0; e < _nGP.size(); ++e) _ofs[e + 1] = _ofs[e] + _nGP[e];
        _data.resize(_ofs.back());
        _final = true;
    }

    // -------- sizes --------
    std::size_t nElements() const { return _nGP.size(); }
    std::size_t nGauss(ElemIdx e) const { assert(e < _nGP.size()); return _nGP[e]; }
    std::size_t totalGP() const { return _data.size(); }

    // -------- access --------
    T&       getElemGP(ElemIdx e, GPIdx g)       { assert(e < _nGP.size() && g < _nGP[e]); return _data[_ofs[e] + g]; }
    const T& getElemGP(ElemIdx e, GPIdx g) const { assert(e < _nGP.size() && g < _nGP[e]); return _data[_ofs[e] + g]; }
    void     setElemGP(ElemIdx e, GPIdx g, const T& v){ assert(e < _nGP.size() && g < _nGP[e]); _data[_ofs[e] + g] = v; }

    // -------- bulk --------
    T*       raw()       { return _data.data(); }
    const T* raw() const { return _data.data(); }

private:
    void _build(const std::vector<std::size_t>& v){
        _nGP = v;
        _ofs.assign(_nGP.size() + 1, 0);
        for (std::size_t e = 0; e < _nGP.size(); ++e) _ofs[e + 1] = _ofs[e] + _nGP[e];
        _data.resize(_ofs.back());
        _final = true;
    }

    std::vector<std::size_t> _nGP;   // per-element GP counts
    std::vector<std::size_t> _ofs;   // prefix sums (size nElem+1)
    std::vector<T>           _data;  // flat storage (size sum_e nGP[e])
    bool _final = false;             // finalized shape flag
};


// =========================== FRAMES: REAL TENSORS ============================
// DeformationFrame
// Purpose: Per-time deformation gradients F(e,g).
// StressFrame
// Purpose: Per-time stresses σ(e,g) and P(e,g) (same ragged layout).
// ============================================================================
struct DeformationFrame {
    RaggedElemField<mat3d> F;
    void setElemShape(const std::vector<std::size_t>& nGPperElem){ F.setElemShape(nGPperElem); }
    // Two-pass: F.prepare(nElem); F.setGaussCount(e,ng); F.finalize();
};

struct StressFrame {
    RaggedElemField<mat3d> sigma;
    RaggedElemField<mat3d> P;
    void setElemShape(const std::vector<std::size_t>& nGPperElem){
        sigma.setElemShape(nGPperElem);
        P    .setElemShape(nGPperElem);
    }
    // Two-pass available via sigma/P.prepare(...), setGaussCount(...), finalize()
};

// ===================== SERIES: REAL (MEASURED/COMPUTED) =====================
// TensorSeries<TFrame>
// Purpose:
//   Time series of one tensor field kind with container-level shape cache.
//   New frames inherit the cached element×GP shape.
// Compatible with DeformationFrame and StressFrame via specializations below.
// ============================================================================

// Generic one-channel series (e.g., deformations F)
class Deformations {
public:
    // ---- shape cache for all frames ----
    void setElemShape(const std::vector<std::size_t>& nGPperElem){ _nGPperElem = nGPperElem; _elemReady = true; _applyShapeAll(); }
    void beginElemShape(std::size_t nElem){ _nGPperElem.assign(nElem, 0); _elemReady = false; }
    void setElemGaussCount(ElemIdx e, std::size_t ng){ assert(e < _nGPperElem.size()); _nGPperElem[e] = ng; }
    void finalizeElemShape(){ _elemReady = true; _applyShapeAll(); }

    // ---- time management ----
    TimeIdx addTime(){
        auto t = series.addTime();
        if (_elemReady) series.getTime(t).setElemShape(_nGPperElem);
        return t;
    }

    // ---- data access ----
    void   setF(TimeIdx t, ElemIdx e, GPIdx g, const mat3d& M){ series.getTime(t).F.setElemGP(e, g, M); }
    mat3d& refF(TimeIdx t, ElemIdx e, GPIdx g)               { return series.getTime(t).F.getElemGP(e, g); }
    const mat3d& crefF(TimeIdx t, ElemIdx e, GPIdx g) const  { return series.getTime(t).F.getElemGP(e, g); }

    // ---- sizes ----
    std::size_t nTimes() const { return series.nTimes(); }
    std::size_t nElements(TimeIdx t) const { return series.getTime(t).F.nElements(); }
    std::size_t nGauss(TimeIdx t, ElemIdx e) const { return series.getTime(t).F.nGauss(e); }
    std::size_t totalGP(TimeIdx t) const { return series.getTime(t).F.totalGP(); }

    TimeSeries<DeformationFrame> series;

private:
    void _applyShapeAll(){
        for (std::size_t k = 0; k < series.nTimes(); ++k)
            series.getTime((TimeIdx)k).setElemShape(_nGPperElem);
    }

    std::vector<std::size_t> _nGPperElem;
    bool _elemReady = false;
};

// Two-channel series (σ and P) with shared shape
class Stresses {
public:
    // ---- shape cache ----
    void setElemShape(const std::vector<std::size_t>& nGPperElem){ _nGPperElem = nGPperElem; _elemReady = true; _applyShapeAll(); }
    void beginElemShape(std::size_t nElem){ _nGPperElem.assign(nElem, 0); _elemReady = false; }
    void setElemGaussCount(ElemIdx e, std::size_t ng){ assert(e < _nGPperElem.size()); _nGPperElem[e] = ng; }
    void finalizeElemShape(){ _elemReady = true; _applyShapeAll(); }

    // ---- time ----
    TimeIdx addTime(){
        auto t = series.addTime();
        if (_elemReady) series.getTime(t).setElemShape(_nGPperElem);
        return t;
    }

    // ---- access ----
    void   setSigma(TimeIdx t, ElemIdx e, GPIdx g, const mat3d& S){ series.getTime(t).sigma.setElemGP(e, g, S); }
    void   setP    (TimeIdx t, ElemIdx e, GPIdx g, const mat3d& Pm){ series.getTime(t).P    .setElemGP(e, g, Pm); }
    mat3d& refSigma(TimeIdx t, ElemIdx e, GPIdx g)                { return series.getTime(t).sigma.getElemGP(e, g); }
    mat3d& refP    (TimeIdx t, ElemIdx e, GPIdx g)                { return series.getTime(t).P    .getElemGP(e, g); }
    const mat3d& crefSigma(TimeIdx t, ElemIdx e, GPIdx g) const   { return series.getTime(t).sigma.getElemGP(e, g); }
    const mat3d& crefP    (TimeIdx t, ElemIdx e, GPIdx g) const   { return series.getTime(t).P    .getElemGP(e, g); }

    // ---- sizes ----
    std::size_t nTimes() const { return series.nTimes(); }
    std::size_t nElements(TimeIdx t) const { return series.getTime(t).sigma.nElements(); }
    std::size_t nGauss(TimeIdx t, ElemIdx e) const { return series.getTime(t).sigma.nGauss(e); }
    std::size_t totalGP(TimeIdx t) const { return series.getTime(t).sigma.totalGP(); }

    TimeSeries<StressFrame> series;

private:
    void _applyShapeAll(){
        for (std::size_t k = 0; k < series.nTimes(); ++k)
            series.getTime((TimeIdx)k).setElemShape(_nGPperElem);
    }

    std::vector<std::size_t> _nGPperElem;
    bool _elemReady = false;
};

// ====================== VF-COLLECTION: VIRTUAL TENSORS =======================
// VirtualDeformationFrame
// Purpose: Same as DeformationFrame but scoped to a single virtual field.
// VirtualDeformations
// Purpose: Collection of nVF time series, each with the SAME element×GP shape.
// API mirrors the displacement-side VirtualFields: VF-major → time.
// ============================================================================
struct VirtualDeformationFrame {
    RaggedElemField<mat3d> F;
    void setElemShape(const std::vector<std::size_t>& nGPperElem){ F.setElemShape(nGPperElem); }
};

class VirtualDeformations {
public:
    // ---- virtual-field management ----
    void resizeVF(std::size_t nVF){ _vf.resize(nVF); }
    std::size_t nVF() const { return _vf.size(); }

    // ---- shared shape cache across all VFs ----
    void setElemShape(const std::vector<std::size_t>& nGPperElem){
        _nGPperElem = nGPperElem;
        _elemReady = true;
        _applyShapeAll();
    }
    void beginElemShape(std::size_t nElem){ _nGPperElem.assign(nElem, 0); _elemReady = false; }
    void setElemGaussCount(ElemIdx e, std::size_t ng){ assert(e < _nGPperElem.size()); _nGPperElem[e] = ng; }
    void finalizeElemShape(){ _elemReady = true; _applyShapeAll(); }

    // ---- time per VF ----
    TimeIdx addTime(VFIdx v){
        auto& ts = _vf[(std::size_t)v];
        auto t = ts.addTime();
        if (_elemReady) ts.getTime(t).setElemShape(_nGPperElem);
        return t;
    }

    // ---- access ----
    void   setF(VFIdx v, TimeIdx t, ElemIdx e, GPIdx g, const mat3d& M){ _vf[(std::size_t)v].getTime(t).F.setElemGP(e, g, M); }
    mat3d& refF(VFIdx v, TimeIdx t, ElemIdx e, GPIdx g)               { return _vf[(std::size_t)v].getTime(t).F.getElemGP(e, g); }
    const mat3d& crefF(VFIdx v, TimeIdx t, ElemIdx e, GPIdx g) const  { return _vf[(std::size_t)v].getTime(t).F.getElemGP(e, g); }

    // ---- sizes ----
    std::size_t nTimes(VFIdx v) const { return _vf[(std::size_t)v].nTimes(); }
    std::size_t nElements(VFIdx v, TimeIdx t) const { return _vf[(std::size_t)v].getTime(t).F.nElements(); }
    std::size_t nGauss   (VFIdx v, TimeIdx t, ElemIdx e) const { return _vf[(std::size_t)v].getTime(t).F.nGauss(e); }
    std::size_t totalGP  (VFIdx v, TimeIdx t) const { return _vf[(std::size_t)v].getTime(t).F.totalGP(); }

    // optional direct series access
    TimeSeries<VirtualDeformationFrame>&       getVF(VFIdx v)       { return _vf[(std::size_t)v]; }
    const TimeSeries<VirtualDeformationFrame>& getVF(VFIdx v) const { return _vf[(std::size_t)v]; }

private:
    void _applyShapeAll(){
        for (auto& ts : _vf)
            for (std::size_t k = 0; k < ts.nTimes(); ++k)
                ts.getTime((TimeIdx)k).setElemShape(_nGPperElem);
    }

    std::vector<std::size_t> _nGPperElem; // shared elem→GP shape across all VFs
    bool _elemReady = false;

    std::vector< TimeSeries<VirtualDeformationFrame> > _vf; // VF-major → time
};
