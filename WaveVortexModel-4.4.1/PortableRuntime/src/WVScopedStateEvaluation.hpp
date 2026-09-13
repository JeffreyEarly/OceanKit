#pragma once

#include "WaveVortexKernel/WVKernelTypes.hpp"

namespace wavevortex::runtime::detail {

// Join only an explicitly active, immutable evaluation of this exact state.
// Otherwise own a fresh scope and close it on every return or exception.
template<class Engine>
class WVScopedStateEvaluation final {
public:
  WVScopedStateEvaluation(Engine& engine,const WVState& state):engine_(engine) {
    if(engine_.stateEvaluationActive()) status_=engine_.validateStateEvaluation(state);
    else {status_=engine_.beginStateEvaluation(state); owns_=static_cast<bool>(status_);}
  }
  ~WVScopedStateEvaluation() {if(owns_) engine_.endStateEvaluation();}
  WVScopedStateEvaluation(const WVScopedStateEvaluation&)=delete;
  WVScopedStateEvaluation& operator=(const WVScopedStateEvaluation&)=delete;
  const WVKernelStatus& status() const noexcept {return status_;}
private:
  Engine& engine_;
  bool owns_=false;
  WVKernelStatus status_=WVKernelStatus::ok();
};

} // namespace wavevortex::runtime::detail
