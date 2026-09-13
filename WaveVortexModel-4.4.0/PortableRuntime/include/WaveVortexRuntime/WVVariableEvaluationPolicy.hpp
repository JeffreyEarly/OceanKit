#pragma once
#include <cstdint>

namespace wavevortex::runtime {
enum class WVVariableEvaluationPolicy : std::uint8_t { reuse, lowMemory };

inline const char* variableEvaluationPolicyIdentifier(WVVariableEvaluationPolicy policy) noexcept {
  return policy == WVVariableEvaluationPolicy::lowMemory ? "low-memory" : "reuse";
}

} // namespace wavevortex::runtime
