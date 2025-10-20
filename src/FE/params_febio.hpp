// fe/params_febio.hpp
#pragma once
#include <vector>
#include <string>
#include <FECore/FEModel.h>
#include <FECore/FEParam.h>
#include "state/vfm_state.hpp"

class FEBioParameterApplier final {
public:
  FEBioParameterApplier(FEModel& fem, VFMState& s): _fem(fem), _s(s) { _build_cache(); }

  // p.size() must equal _s.params.size()
  bool apply(const std::vector<double>& p, std::string& err) {
    if (p.size()!=_s.params.size()){ err="param size mismatch"; return false; }
    for (size_t i=0;i<p.size();++i){
      if (_ptrs[i]==nullptr){ err="null FE ptr for "+_s.params[i].spec.name; return false; }
      *(_ptrs[i]) = p[i];          // write into FEBio
      _s.params[i].value = p[i];   // mirror in state
    }
    return true;
  }

private:
  void _build_cache(){
    _ptrs.assign(_s.params.size(), nullptr);
    for (size_t i=0;i<_s.params.size();++i){
      FEParamValue v = _fem.GetParameterValue(ParamString(_s.params[i].spec.name.c_str()));
      if (v.isValid() && v.type()==FE_PARAM_DOUBLE) _ptrs[i] = static_cast<double*>(v.data_ptr());
    }
  }
  FEModel& _fem; VFMState& _s;
  std::vector<double*> _ptrs;
};
