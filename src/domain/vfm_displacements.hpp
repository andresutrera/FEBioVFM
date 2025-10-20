// vfm_displacements.hpp
#pragma once
#include <vector>
#include <cstddef>
#include <cassert>
#include <FECore/vec3d.h>
#include "domain/vfm_core_series.hpp"

// --- Index aliases kept small and POD-only ---
using NodeIdx = std::size_t; // node index in [0..nNodes)
using TimeIdx = int;         // time/frame index in [0..nTimes)
using VFIdx   = int;         // virtual-field index in [0..nVF)

// ============================================================================
// NodalField<T>
// Purpose:
//   Dense per-node storage for a value type T (e.g., vec3d).
// Guarantees:
//   O(1) get/set by NodeIdx; size fixed after resizeNodes() until next resize.
// Notes:
//   No implicit allocation on access; caller controls size explicitly.
// ============================================================================
template <typename T>
class NodalField {
public:
    NodalField() = default;

    // ---- shape / size ----
    void         resizeNodes(std::size_t n) { _d.resize(n); }
    std::size_t  size() const               { return _d.size(); }

    // ---- accessors ----
    T&           getNode(NodeIdx i)         { assert(i < _d.size()); return _d[i]; }
    const T&     getNode(NodeIdx i) const   { assert(i < _d.size()); return _d[i]; }
    void         setNode(NodeIdx i, const T& v){ assert(i < _d.size()); _d[i] = v; }

    // ---- bulk ----
    T*           raw()                      { return _d.data(); }
    const T*     raw() const                { return _d.data(); }

private:
    std::vector<T> _d;
};

// ============================================================================
// Frames (single time)
// MeasuredFrame  : experimental nodal displacements u_i(t)
// VirtualFrame   : virtual-field nodal displacements u^v_i(t)
// LoadFrame      : measured nodal forces F_i(t)
// ============================================================================
struct MeasuredFrame {
    NodalField<vec3d> u;
    void setNodalSize(std::size_t n) { u.resizeNodes(n); }
};

struct VirtualFrame {
    NodalField<vec3d> u;
    void setNodalSize(std::size_t n) { u.resizeNodes(n); }
};

struct LoadFrame {
    NodalField<vec3d> F;
    void setNodalSize(std::size_t n) { F.resizeNodes(n); }
};

// ============================================================================
// MeasuredData
// Purpose:
//   Time series of measured nodal displacements u_i(t) with container-level
//   node-count cache. New frames inherit the cached node size.
// API:
//   setNodalSize(n), addTime(), setU/refU/crefU.
// ============================================================================
class MeasuredData {
public:
    // ---- shape cache for all frames ----
    void setNodalSize(std::size_t nNodes) { _nNodes=nNodes; _nodalReady=true; _applyNodalToAll(); }

    // ---- time management ----
    TimeIdx addTime() {
        auto t = series.addTime();
        if (_nodalReady) series.getTime(t).setNodalSize(_nNodes);
        return t;
    }

    // ---- accessors ----
    void        setU (TimeIdx t, NodeIdx i, const vec3d& v) { series.getTime(t).u.setNode(i, v); }
    vec3d&      refU (TimeIdx t, NodeIdx i)                 { return series.getTime(t).u.getNode(i); }
    const vec3d&crefU(TimeIdx t, NodeIdx i) const           { return series.getTime(t).u.getNode(i); }

    // ---- sizes ----
    std::size_t nTimes() const { return series.nTimes(); }

    TimeSeries<MeasuredFrame> series; // exposed for iteration when needed

private:
    void _applyNodalToAll() {
        for (std::size_t k=0; k<series.nTimes(); ++k) series.getTime((TimeIdx)k).setNodalSize(_nNodes);
    }
    std::size_t _nNodes=0;
    bool _nodalReady=false;
};

// ============================================================================
// MeasuredLoad
// Purpose:
//   Time series of measured nodal loads F_i(t) with shared node size cache.
// API mirrors MeasuredData for symmetry.
// ============================================================================
class MeasuredLoad {
public:
    void setNodalSize(std::size_t nNodes) { _nNodes=nNodes; _nodalReady=true; _applyNodalToAll(); }

    TimeIdx addTime() {
        auto t = series.addTime();
        if (_nodalReady) series.getTime(t).setNodalSize(_nNodes);
        return t;
    }

    void        setF (TimeIdx t, NodeIdx i, const vec3d& v) { series.getTime(t).F.setNode(i, v); }
    vec3d&      refF (TimeIdx t, NodeIdx i)                 { return series.getTime(t).F.getNode(i); }
    const vec3d&crefF(TimeIdx t, NodeIdx i) const           { return series.getTime(t).F.getNode(i); }

    std::size_t nTimes() const { return series.nTimes(); }

    TimeSeries<LoadFrame> series;

private:
    void _applyNodalToAll() {
        for (std::size_t k=0; k<series.nTimes(); ++k) series.getTime((TimeIdx)k).setNodalSize(_nNodes);
    }
    std::size_t _nNodes=0;
    bool _nodalReady=false;
};

// ============================================================================
// VirtualFields (displacements)
// Purpose:
//   Collection of nVF virtual fields, each a time series of nodal vectors,
//   sharing a single node-count cache. VF-major → time indexing.
// API:
//   resizeVF(nVF), setNodalSize(n), addTime(v), setU/refU/crefU.
// ============================================================================
class VirtualFields {
public:
    // ---- VF management ----
    void        resizeVF(std::size_t nVF) { _vf.resize(nVF); }
    std::size_t nVF() const               { return _vf.size(); }

    // ---- shape cache shared by all VFs ----
    void setNodalSize(std::size_t nNodes) { _nNodes=nNodes; _nodalReady=true; _applyNodalToAll(); }

    // ---- time per VF ----
    TimeIdx addTime(VFIdx v) {
        auto& ts = _vf[(std::size_t)v];
        auto t = ts.addTime();
        if (_nodalReady) ts.getTime(t).setNodalSize(_nNodes);
        return t;
    }

    // ---- accessors ----
    void        setU (VFIdx v, TimeIdx t, NodeIdx i, const vec3d& val) { _vf[(std::size_t)v].getTime(t).u.setNode(i, val); }
    vec3d&      refU (VFIdx v, TimeIdx t, NodeIdx i)                   { return _vf[(std::size_t)v].getTime(t).u.getNode(i); }
    const vec3d&crefU(VFIdx v, TimeIdx t, NodeIdx i) const             { return _vf[(std::size_t)v].getTime(t).u.getNode(i); }

    // ---- direct access to a VF series (optional) ----
    TimeSeries<VirtualFrame>&       getVF(VFIdx v)       { return _vf[(std::size_t)v]; }
    const TimeSeries<VirtualFrame>& getVF(VFIdx v) const { return _vf[(std::size_t)v]; }

private:
    void _applyNodalToAll() {
        for (auto& ts : _vf)
            for (std::size_t k=0; k<ts.nTimes(); ++k)
                ts.getTime((TimeIdx)k).setNodalSize(_nNodes);
    }

    std::size_t _nNodes=0;
    bool _nodalReady=false;

    std::vector< TimeSeries<VirtualFrame> > _vf; // VF-major → time
};
