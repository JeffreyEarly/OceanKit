#pragma once

#include "WaveVortexRuntime/WVModel.hpp"

namespace wavevortex::runtime::detail {

// Construction-resolved transform boundary owned by WVModel. The generic
// model never retains a transform-specific side pointer or performs a
// downcast; named adapters own their concrete numerical systems.
class WVResolvedModelSystem {
public:
  virtual ~WVResolvedModelSystem() = default;
  virtual WVIntegrationSystem &integrationSystem() noexcept = 0;
  virtual const WVIntegrationSystem &integrationSystem() const noexcept = 0;
  virtual const std::string &forcingScheduleIdentifier() const noexcept = 0;
  virtual const std::string &kernelProviderIdentifier() const noexcept = 0;
  virtual const std::string &kernelProviderLibraryIdentity() const noexcept = 0;
  virtual WVKernelStatus
  initializeObserverState(WVMutableIntegrationState &state) const = 0;
  virtual void populateMetrics(WVModelMetrics &metrics) const noexcept = 0;
  virtual WVTransformConstantStratificationKernel *
  constantStratificationKernel() noexcept {
    return nullptr;
  }
};

WVKernelStatus createConstantStratificationModelSystem(
    const WVTransformConstantStratificationConfiguration &configuration,
    const WVFrozenForcingSchedule &schedule,
    const WVPortableObserverDescriptor *descriptor,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVResolvedModelSystem> &system);

WVKernelStatus createBarotropicQGModelSystem(
    const WVTransformBarotropicQGConfiguration &configuration,
    const WVFrozenForcingSchedule &schedule,
    const WVPortableObserverDescriptor *descriptor,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVResolvedModelSystem> &system);

WVKernelStatus createPersistedModelSystem(
    const WVCheckpointInspection &inspection,
    const WVFrozenForcingSchedule &schedule,
    const WVPortableObserverDescriptor *descriptor,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    std::unique_ptr<WVFFTEngine> engine,
    std::unique_ptr<WVResolvedModelSystem> &system);
WVKernelStatus validatePersistedModelForcingSchedule(
    const WVCheckpointInspection &inspection,
    const WVFrozenForcingSchedule &schedule,
    const WVExtensionCatalog &catalog);
WVCheckpointInspection modelCheckpointInspection(
    const WVCheckpoint &checkpoint);
const WVTransformConstantStratificationConfiguration *
legacyModelOutputPlanningConfiguration(
    const WVCheckpointInspection &inspection) noexcept;

// Checkpoint storage is intentionally resolved through the named transform
// adapter so WVModelState only consumes ordered coefficient-family views.
WVKernelStatus validateModelCheckpointState(
    WVCheckpoint &checkpoint, const WVIntegrationStateLayout &layout);
WVComplex64 *modelCheckpointCoefficientData(
    WVCheckpoint &checkpoint, const WVIntegrationStateLayout &layout,
    std::size_t family) noexcept;
const WVComplex64 *modelCheckpointCoefficientData(
    const WVCheckpoint &checkpoint, const WVIntegrationStateLayout &layout,
    std::size_t family) noexcept;
WVMutableState modelCheckpointLegacyView(WVCheckpoint &checkpoint) noexcept;
void setModelCheckpointTimes(WVCheckpoint &checkpoint, double t,
                             double t0) noexcept;

} // namespace wavevortex::runtime::detail
