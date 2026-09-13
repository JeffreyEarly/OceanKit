#pragma once

#include "WaveVortexKernel/WVKernelTypes.hpp"
#include "WaveVortexRuntime/WVVariableEvaluationPolicy.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <new>
#include <utility>
#include <vector>

namespace wavevortex::runtime {

enum class WVVariableEvaluationNode : std::uint8_t {
  registeredVariable, phaseFactors, physicalField, derivative, reduction,
  componentCoefficients, projection, forcingTendency, gridCalculus
};

// Owner and immutable state are bound by the enclosing evaluation. Geometry is
// zero for natural-grid quantities; only sampling operations use a geometry key.
struct WVVariableEvaluationKey {
  WVVariableEvaluationNode node = WVVariableEvaluationNode::registeredVariable;
  std::uint32_t variable = 0, component = 0, derivative = 0, reference = 0;
  std::uint32_t forcingOrdinal = 0, stage = 0;
  std::uint64_t geometry = 0;
  bool operator==(const WVVariableEvaluationKey& b) const noexcept {
    return node==b.node && variable==b.variable && component==b.component &&
        derivative==b.derivative && reference==b.reference &&
        forcingOrdinal==b.forcingOrdinal && stage==b.stage && geometry==b.geometry;
  }
};

struct WVVariableEvaluationMetrics {
  std::size_t contexts = 0, producerExecutions = 0, cacheHits = 0;
  std::size_t evictions = 0, recomputations = 0, duplicateExecutions = 0;
  std::size_t liveBytes = 0, highWaterBytes = 0;
};

// Counters taken at numerical producers, independently of the execution ledger.
// Field indices follow each transform's field enum; derivative zero is value.
struct WVVariableProducerMetrics {
  std::size_t stateValidations = 0, phasePreparations = 0;
  std::size_t derivedValidations = 0;
  std::size_t coefficientAssemblies = 0, verticalPreparations = 0;
  std::size_t verticalOperatorExecutions = 0;
  std::size_t verticalMatrixGroupExecutions = 0;
  std::size_t horizontalSpectrumReuses = 0, preparedVerticalDerivatives = 0;
  std::size_t derivativeAdvectionConsumers = 0;
  std::size_t tiledNonlinearExecutions = 0, tiledColumnInverses = 0;
  std::size_t tiledRowInverses = 0, tiledReusedColumns = 0;
  std::size_t horizontalSpeedReductions = 0;
  std::size_t verticalSpeedReductions = 0;
  std::size_t energyReductions = 0;
  std::array<std::size_t,4> tendencyReconstructions{};
  std::array<std::array<std::array<std::size_t,5>,4>,16> reconstructions{};
};

// Execution ledger for buffers owned by a family adapter or output workspace.
// All callers use the same dependency keys and lifecycle. This ledger never
// owns scientific state, decides that a mutable pointer is unchanged, or moves
// an adapter's arrays. Prepare keys before entering a time-stepping loop.
class WVVariableEvaluationContext final {
public:
  WVKernelStatus prepare(const std::vector<WVVariableEvaluationKey>& keys) {
    if(active_) return invalid("Cannot change an active variable evaluation plan.");
    try {
      entries_.reserve(entries_.size()+keys.size());
      for(const auto& key:keys) if(!find(key)) entries_.push_back({key});
    } catch(const std::bad_alloc&) {
      return {WVKernelStatusCode::allocationFailure,"Unable to prepare variable evaluation entries."};
    }
    return WVKernelStatus::ok();
  }
  WVKernelStatus begin(const void* owner, WVVariableEvaluationPolicy policy=WVVariableEvaluationPolicy::reuse) {
    if(active_) return {WVKernelStatusCode::reentrantExecution,"Variable evaluation is already active."};
    if(!owner) return invalid("Variable evaluation requires an owner.");
    if(policy!=WVVariableEvaluationPolicy::reuse && policy!=WVVariableEvaluationPolicy::lowMemory)
      return invalid("Unknown variable evaluation policy.");
    for(auto& entry:entries_) {entry.state=State::empty; entry.executions=0; entry.bytes=0; entry.pins=0;}
    owner_=owner; policy_=policy; active_=true; ++generation_; ++metrics_.contexts;
    metrics_.liveBytes=0;
    return WVKernelStatus::ok();
  }
  void end() noexcept {
    for(auto& entry:entries_) {entry.state=State::empty; entry.pins=0; entry.bytes=0;}
    metrics_.liveBytes=0; active_=false; owner_=nullptr;
  }
  bool active() const noexcept {return active_;}
  std::uint64_t generation() const noexcept {return generation_;}
  WVVariableEvaluationPolicy policy() const noexcept {return policy_;}
  const WVVariableEvaluationMetrics& metrics() const noexcept {return metrics_;}
  std::size_t persistentBytes() const noexcept {return entries_.capacity()*sizeof(Entry);}
  bool ready(const WVVariableEvaluationKey& key) const noexcept {
    const auto* entry=find(key); return active_ && entry && entry->state==State::ready;
  }
  WVKernelStatus validate(const void* owner, std::uint64_t generation) const {
    return active_ && owner==owner_ && generation==generation_ ? WVKernelStatus::ok() :
        invalid("Variable evaluation view is stale or belongs to another owner.");
  }
  template<class Producer>
  WVKernelStatus evaluate(const WVVariableEvaluationKey& key, std::size_t bytes, Producer&& producer) {
    if(!active_) return invalid("Variable evaluation requires an active immutable-state scope.");
    auto* entry=find(key);
    if(!entry) return invalid("Variable evaluation node was not prepared.");
    if(entry->state==State::ready) {
      if(entry->bytes!=bytes) return invalid("Cached variable extent differs from the prepared request.");
      ++metrics_.cacheHits; return WVKernelStatus::ok();
    }
    if(entry->state==State::computing) return invalid("Cyclic variable evaluation dependency.");
    entry->state=State::computing;
    // Restore an empty node after failures, including exceptions from a producer.
    struct Guard {Entry& entry; ~Guard(){if(entry.state==State::computing) entry.state=State::empty;}} guard{*entry};
    const auto status=producer();
    if(!status) return status;
    if(entry->executions) {
      ++metrics_.recomputations;
      if(policy_==WVVariableEvaluationPolicy::reuse) ++metrics_.duplicateExecutions;
    }
    ++entry->executions; ++metrics_.producerExecutions;
    entry->state=State::ready; entry->bytes=bytes;
    metrics_.liveBytes+=bytes;
    metrics_.highWaterBytes=std::max(metrics_.highWaterBytes,metrics_.liveBytes);
    return WVKernelStatus::ok();
  }
  // One fused producer can materialize several prepared nodes atomically. A
  // group is either entirely cached or entirely empty; callers using the
  // low-memory policy must evict every ready member before recomputing it.
  template<class Producer>
  WVKernelStatus evaluateGroup(
      const std::vector<std::pair<WVVariableEvaluationKey,std::size_t>>& nodes,
      Producer&& producer) {
    if(!active_) return invalid("Variable evaluation requires an active immutable-state scope.");
    if(nodes.empty()) return invalid("A fused variable evaluation group cannot be empty.");
    std::size_t readyCount=0;
    for(std::size_t index=0;index<nodes.size();++index) {
      const auto& node=nodes[index];
      auto* entry=find(node.first);
      if(!entry) return invalid("Fused variable evaluation node was not prepared.");
      for(std::size_t previous=0;previous<index;++previous)
        if(nodes[previous].first==node.first)
          return invalid("Fused variable evaluation contains a duplicate node.");
      if(entry->state==State::computing)
        return invalid("Cyclic fused variable evaluation dependency.");
      if(entry->state==State::ready) {
        if(entry->bytes!=node.second)
          return invalid("Cached fused variable extent differs from the prepared request.");
        ++readyCount;
      }
    }
    if(readyCount==nodes.size()) {
      metrics_.cacheHits+=nodes.size();
      return WVKernelStatus::ok();
    }
    if(readyCount)
      return invalid("Fused variable evaluation overlaps already-ready nodes; evict the complete group before recomputation.");
    for(const auto& node:nodes) find(node.first)->state=State::computing;
    struct Guard {
      WVVariableEvaluationContext& context;
      const std::vector<std::pair<WVVariableEvaluationKey,std::size_t>>& nodes;
      ~Guard(){for(const auto& node:nodes) {
        auto* entry=context.find(node.first);
        if(entry->state==State::computing) entry->state=State::empty;
      }}
    } guard{*this,nodes};
    const auto status=producer();
    if(!status) return status;
    bool repeated=false;
    for(const auto& node:nodes) repeated=repeated || find(node.first)->executions!=0;
    if(repeated) {
      ++metrics_.recomputations;
      if(policy_==WVVariableEvaluationPolicy::reuse) ++metrics_.duplicateExecutions;
    }
    ++metrics_.producerExecutions;
    for(const auto& node:nodes) {
      auto* entry=find(node.first);
      ++entry->executions;
      entry->state=State::ready;
      entry->bytes=node.second;
      metrics_.liveBytes+=entry->bytes;
    }
    metrics_.highWaterBytes=std::max(metrics_.highWaterBytes,metrics_.liveBytes);
    return WVKernelStatus::ok();
  }
  bool pin(const WVVariableEvaluationKey& key) noexcept {
    auto* entry=find(key); if(!active_ || !entry || entry->state!=State::ready) return false;
    ++entry->pins; return true;
  }
  void unpin(const WVVariableEvaluationKey& key) noexcept {
    if(auto* entry=find(key); entry && entry->pins) --entry->pins;
  }
  // Call only at a proven last-consumer boundary. The adapter can release or
  // reuse the underlying storage only after this returns true.
  bool evict(const WVVariableEvaluationKey& key) noexcept {
    auto* entry=find(key);
    if(!active_ || policy_!=WVVariableEvaluationPolicy::lowMemory || !entry ||
        entry->state!=State::ready || entry->pins) return false;
    metrics_.liveBytes-=entry->bytes; entry->bytes=0; entry->state=State::empty;
    ++metrics_.evictions; return true;
  }
private:
  enum class State : std::uint8_t {empty, computing, ready};
  struct Entry {WVVariableEvaluationKey key; State state=State::empty; std::size_t executions=0,bytes=0,pins=0;};
  static WVKernelStatus invalid(const char* message) {return {WVKernelStatusCode::invalidConfiguration,message};}
  Entry* find(const WVVariableEvaluationKey& key) noexcept {
    for(auto& entry:entries_) {
      if(entry.key==key) return &entry;
    }
    return nullptr;
  }
  const Entry* find(const WVVariableEvaluationKey& key) const noexcept {
    for(const auto& entry:entries_) {
      if(entry.key==key) return &entry;
    }
    return nullptr;
  }
  std::vector<Entry> entries_;
  const void* owner_=nullptr;
  std::uint64_t generation_=0;
  bool active_=false;
  WVVariableEvaluationPolicy policy_=WVVariableEvaluationPolicy::reuse;
  WVVariableEvaluationMetrics metrics_;
};

} // namespace wavevortex::runtime
