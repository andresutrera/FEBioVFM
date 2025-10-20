// vfm.hpp
#pragma once
#include <vector>
#include <cstddef>
#include <cassert>
#include <cstring> // optional for bulk memcpy
#include <FECore/mat3d.h>

using NodeIdx = std::size_t;
using ElemIdx = std::size_t;
using GPIdx = std::size_t;
using TimeIdx = int;
using VFIdx = int;

// -------- fields --------
template <typename T>
class NodalField
{
public:
	NodalField() = default;
	void resizeNodes(std::size_t n) { _d.resize(n); }
	std::size_t size() const { return _d.size(); }

	// access
	T &getNode(NodeIdx i)
	{
		assert(i < _d.size());
		return _d[i];
	}
	const T &getNode(NodeIdx i) const
	{
		assert(i < _d.size());
		return _d[i];
	}
	void setNode(NodeIdx i, const T &v)
	{
		assert(i < _d.size());
		_d[i] = v;
	}

	// bulk
	T *raw() { return _d.data(); }
	const T *raw() const { return _d.data(); }

private:
	std::vector<T> _d;
};

template <typename T>
class RaggedElemField
{
public:
	RaggedElemField() = default;

	// one-shot shape
	void setElemShape(const std::vector<std::size_t> &nGPperElem) { _buildFromVector(nGPperElem); }

	// two-pass incremental shape
	void prepare(std::size_t nElem)
	{
		_nGP.assign(nElem, 0);
		_ofs.assign(nElem + 1, 0);
		_data.clear();
		_finalized = false;
	}
	void setGaussCount(ElemIdx e, std::size_t ng)
	{
		assert(!_finalized && e < _nGP.size());
		_nGP[e] = ng;
	}
	void finalize()
	{
		assert(!_finalized);
		for (std::size_t e = 0; e < _nGP.size(); ++e)
			_ofs[e + 1] = _ofs[e] + _nGP[e];
		_data.resize(_ofs.back());
		_finalized = true;
	}

	// sizes
	std::size_t nElements() const { return _nGP.size(); }
	std::size_t nGauss(ElemIdx e) const
	{
		assert(e < _nGP.size());
		return _nGP[e];
	}
	std::size_t totalGP() const { return _data.size(); }

	// access
	T &getElemGP(ElemIdx e, GPIdx g)
	{
		assert(e < _nGP.size() && g < _nGP[e]);
		return _data[_ofs[e] + g];
	}
	const T &getElemGP(ElemIdx e, GPIdx g) const
	{
		assert(e < _nGP.size() && g < _nGP[e]);
		return _data[_ofs[e] + g];
	}
	void setElemGP(ElemIdx e, GPIdx g, const T &v)
	{
		assert(e < _nGP.size() && g < _nGP[e]);
		_data[_ofs[e] + g] = v;
	}

	// bulk
	T *raw() { return _data.data(); }
	const T *raw() const { return _data.data(); }

private:
	void _buildFromVector(const std::vector<std::size_t> &v)
	{
		_nGP = v;
		_ofs.assign(_nGP.size() + 1, 0);
		for (std::size_t e = 0; e < _nGP.size(); ++e)
			_ofs[e + 1] = _ofs[e] + _nGP[e];
		_data.resize(_ofs.back());
		_finalized = true;
	}
	std::vector<std::size_t> _nGP, _ofs;
	std::vector<T> _data;
	bool _finalized = false;
};

// -------- frames (single time) --------
struct MeasuredFrame
{
	NodalField<vec3d> u;
	RaggedElemField<mat3d> F;

	// shape setters (independent)
	void setNodalSize(std::size_t nNodes) { u.resizeNodes(nNodes); }
	// element shape: either one-shot or two-pass via F.prepare/setGaussCount/finalize
	void setElemShape(const std::vector<std::size_t> &nGPperElem) { F.setElemShape(nGPperElem); }
};

struct VirtualFrame
{
	NodalField<vec3d> u;
	RaggedElemField<mat3d> F;
	void setNodalSize(std::size_t nNodes) { u.resizeNodes(nNodes); }
	void setElemShape(const std::vector<std::size_t> &nGPperElem) { F.setElemShape(nGPperElem); }
};

struct StressFrame
{
	RaggedElemField<mat3d> sigma;
	RaggedElemField<mat3d> P;
	void setElemShape(const std::vector<std::size_t> &nGPperElem)
	{
		sigma.setElemShape(nGPperElem);
		P.setElemShape(nGPperElem);
	}
	// two-pass exposed via sigma/P.prepare(..), setGaussCount(..), finalize()
};

struct LoadFrame
{
	NodalField<vec3d> F;
	void setNodalSize(std::size_t nNodes) { F.resizeNodes(nNodes); }
};

// -------- time series --------
template <typename FrameT>
class TimeSeries
{
public:
	std::size_t nTimes() const { return _t.size(); }
	TimeIdx addTime()
	{
		_t.emplace_back();
		return (TimeIdx)_t.size() - 1;
	}
	FrameT &getTime(TimeIdx t)
	{
		assert(t >= 0 && (std::size_t)t < _t.size());
		return _t[(std::size_t)t];
	}
	const FrameT &getTime(TimeIdx t) const
	{
		assert(t >= 0 && (std::size_t)t < _t.size());
		return _t[(std::size_t)t];
	}

private:
	std::vector<FrameT> _t;
};

// -------- containers --------

// Measured experimental fields
class MeasuredData
{
public:
	// shapes can be set independently and at any time
	void setNodalSize(std::size_t nNodes)
	{
		_nNodes = nNodes;
		_nodalReady = true;
		_applyNodalToAll();
	}
	void setElemShape(const std::vector<std::size_t> &nGPperElem)
	{
		_nGPperElem = nGPperElem;
		_elemReady = true;
		_applyElemToAll();
	}

	// two-pass element shape, container-level helpers
	void beginElemShape(std::size_t nElem)
	{
		_nGPperElem.assign(nElem, 0);
		_elemReady = false;
	}
	void setElemGaussCount(ElemIdx e, std::size_t ng)
	{
		assert(e < _nGPperElem.size());
		_nGPperElem[e] = ng;
	}
	void finalizeElemShape()
	{
		_elemReady = true;
		_applyElemToAll();
	}

	TimeIdx addTime()
	{
		auto t = series.addTime();
		if (_nodalReady)
			series.getTime(t).setNodalSize(_nNodes);
		if (_elemReady)
			series.getTime(t).setElemShape(_nGPperElem);
		return t;
	}

	// node setters
	void setU(TimeIdx t, NodeIdx i, const vec3d &v) { series.getTime(t).u.setNode(i, v); }
	vec3d &refU(TimeIdx t, NodeIdx i) { return series.getTime(t).u.getNode(i); }

	// defgrad setters
	void setF(TimeIdx t, ElemIdx e, GPIdx g, const mat3d &m) { series.getTime(t).F.setElemGP(e, g, m); }
	mat3d &refF(TimeIdx t, ElemIdx e, GPIdx g) { return series.getTime(t).F.getElemGP(e, g); }

	// quick sizes
	std::size_t nTimes() const { return series.nTimes(); }

	TimeSeries<MeasuredFrame> series;

private:
	void _applyNodalToAll()
	{
		for (std::size_t k = 0; k < series.nTimes(); ++k)
			series.getTime((TimeIdx)k).setNodalSize(_nNodes);
	}
	void _applyElemToAll()
	{
		for (std::size_t k = 0; k < series.nTimes(); ++k)
			series.getTime((TimeIdx)k).setElemShape(_nGPperElem);
	}

	std::size_t _nNodes = 0;
	bool _nodalReady = false;
	std::vector<std::size_t> _nGPperElem;
	bool _elemReady = false;
};

// Stresses (Cauchy and 1st Piola)
class Stresses
{
public:
	void setElemShape(const std::vector<std::size_t> &nGPperElem)
	{
		_nGPperElem = nGPperElem;
		_elemReady = true;
		_applyElemToAll();
	}

	// two-pass build
	void beginElemShape(std::size_t nElem)
	{
		_nGPperElem.assign(nElem, 0);
		_elemReady = false;
	}
	void setElemGaussCount(ElemIdx e, std::size_t ng)
	{
		assert(e < _nGPperElem.size());
		_nGPperElem[e] = ng;
	}
	void finalizeElemShape()
	{
		_elemReady = true;
		_applyElemToAll();
	}

	TimeIdx addTime()
	{
		auto t = series.addTime();
		if (_elemReady)
			series.getTime(t).setElemShape(_nGPperElem);
		return t;
	}

	void setSigma(TimeIdx t, ElemIdx e, GPIdx g, const mat3d &s) { series.getTime(t).sigma.setElemGP(e, g, s); }
	void setP(TimeIdx t, ElemIdx e, GPIdx g, const mat3d &p) { series.getTime(t).P.setElemGP(e, g, p); }
	mat3d &refSigma(TimeIdx t, ElemIdx e, GPIdx g) { return series.getTime(t).sigma.getElemGP(e, g); }
	mat3d &refP(TimeIdx t, ElemIdx e, GPIdx g) { return series.getTime(t).P.getElemGP(e, g); }

	TimeSeries<StressFrame> series;

private:
	void _applyElemToAll()
	{
		for (std::size_t k = 0; k < series.nTimes(); ++k)
			series.getTime((TimeIdx)k).setElemShape(_nGPperElem);
	}
	std::vector<std::size_t> _nGPperElem;
	bool _elemReady = false;
};

// Measured nodal loads
class MeasuredLoad
{
public:
	void setNodalSize(std::size_t nNodes)
	{
		_nNodes = nNodes;
		_nodalReady = true;
		_applyNodalToAll();
	}
	TimeIdx addTime()
	{
		auto t = series.addTime();
		if (_nodalReady)
			series.getTime(t).setNodalSize(_nNodes);
		return t;
	}

	void setF(TimeIdx t, NodeIdx i, const vec3d &v) { series.getTime(t).F.setNode(i, v); }
	vec3d &refF(TimeIdx t, NodeIdx i) { return series.getTime(t).F.getNode(i); }

	TimeSeries<LoadFrame> series;

private:
	void _applyNodalToAll()
	{
		for (std::size_t k = 0; k < series.nTimes(); ++k)
			series.getTime((TimeIdx)k).setNodalSize(_nNodes);
	}
	std::size_t _nNodes = 0;
	bool _nodalReady = false;
};

// Virtual fields: VF-major → time
class VirtualFields
{
public:
	void resizeVF(std::size_t nVF) { _vf.resize(nVF); }

	// shapes cached at container, applied to frames
	void setNodalSize(std::size_t nNodes)
	{
		_nNodes = nNodes;
		_nodalReady = true;
		_applyNodalToAll();
	}
	void setElemShape(const std::vector<std::size_t> &nGPperElem)
	{
		_nGPperElem = nGPperElem;
		_elemReady = true;
		_applyElemToAll();
	}

	// two-pass build for element shape
	void beginElemShape(std::size_t nElem)
	{
		_nGPperElem.assign(nElem, 0);
		_elemReady = false;
	}
	void setElemGaussCount(ElemIdx e, std::size_t ng)
	{
		assert(e < _nGPperElem.size());
		_nGPperElem[e] = ng;
	}
	void finalizeElemShape()
	{
		_elemReady = true;
		_applyElemToAll();
	}

	TimeIdx addTime(VFIdx v)
	{
		auto &ts = _vf[(std::size_t)v];
		auto t = ts.addTime();
		if (_nodalReady)
			ts.getTime(t).setNodalSize(_nNodes);
		if (_elemReady)
			ts.getTime(t).setElemShape(_nGPperElem);
		return t;
	}

	// setters
	void setU(VFIdx v, TimeIdx t, NodeIdx i, const vec3d &val) { _vf[(std::size_t)v].getTime(t).u.setNode(i, val); }
	void setF(VFIdx v, TimeIdx t, ElemIdx e, GPIdx g, const mat3d &val) { _vf[(std::size_t)v].getTime(t).F.setElemGP(e, g, val); }

	// refs
	vec3d &refU(VFIdx v, TimeIdx t, NodeIdx i) { return _vf[(std::size_t)v].getTime(t).u.getNode(i); }
	mat3d &refF(VFIdx v, TimeIdx t, ElemIdx e, GPIdx g) { return _vf[(std::size_t)v].getTime(t).F.getElemGP(e, g); }

	TimeSeries<VirtualFrame> &getVF(VFIdx v) { return _vf[(std::size_t)v]; }
	const TimeSeries<VirtualFrame> &getVF(VFIdx v) const { return _vf[(std::size_t)v]; }

private:
	void _applyNodalToAll()
	{
		for (auto &ts : _vf)
			for (std::size_t k = 0; k < ts.nTimes(); ++k)
				ts.getTime((TimeIdx)k).setNodalSize(_nNodes);
	}
	void _applyElemToAll()
	{
		for (auto &ts : _vf)
			for (std::size_t k = 0; k < ts.nTimes(); ++k)
				ts.getTime((TimeIdx)k).setElemShape(_nGPperElem);
	}

	std::size_t _nNodes = 0;
	bool _nodalReady = false;
	std::vector<std::size_t> _nGPperElem;
	bool _elemReady = false;

	std::vector<TimeSeries<VirtualFrame>> _vf;
};
