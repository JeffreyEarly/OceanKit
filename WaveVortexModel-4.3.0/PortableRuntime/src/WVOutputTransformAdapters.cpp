#include "WVOutputTransformAdapters.hpp"

#include "WaveVortexRuntime/WVObserverContracts.hpp"

#include <algorithm>

namespace wavevortex::runtime::detail {
namespace {

WVKernelStatus invalidShape(std::string message) {
  return {WVKernelStatusCode::invalidShape, std::move(message)};
}

} // namespace

WVKernelStatus createLegacyOutputStateLayout(
    const WVPortableObserverDescriptor &descriptor,
    WVIntegrationStateLayout &layout) {
  const auto &blocks = descriptor.record().stateBlocks;
  const auto canonical =
      std::find_if(blocks.begin(), blocks.end(), [](const auto &block) {
        return block.identifier == "Ap";
      });
  if (canonical == blocks.end() || canonical->dimensions.size() != 2)
    return {WVKernelStatusCode::invalidConfiguration,
            "Legacy output planning requires a canonical [Nj,Nkl] Ap block."};
  return WVIntegrationStateLayout::create(
      {canonical->dimensions[0], canonical->dimensions[1]}, descriptor,
      layout);
}

void bindLegacyOutputStateView(
    const WVIntegrationStateLayout &layout,
    WVCoefficientStateStorage &storage,
    WVMutableIntegrationState &state) noexcept {
  if (!layout.hasLegacyCoefficientTriple())
    return;
  const auto shape = layout.coefficientShape();
  auto *families = storage.mutableFamilies();
  state.waveVortex.coefficients = {{families[0].data, shape},
                                   {families[1].data, shape},
                                   {families[2].data, shape}};
}

WVKernelStatus validateOutputCheckpointTemplate(
    const WVIntegrationStateLayout &layout, const WVCheckpoint &checkpoint) {
  if (layout.hasLegacyCoefficientTriple()) {
    const auto shape = layout.coefficientShape();
    if (checkpoint.state.coefficients.shape.rows != shape.rows ||
        checkpoint.state.coefficients.shape.columns != shape.columns)
      return invalidShape(
          "Checkpoint template and output plan coefficient shapes differ.");
    const auto count = shape.elementCount();
    if (checkpoint.state.coefficients.Ap.size() != count ||
        checkpoint.state.coefficients.Am.size() != count ||
        checkpoint.state.coefficients.A0.size() != count)
      return invalidShape("Checkpoint template coefficient storage is "
                          "incomplete.");
    return WVKernelStatus::ok();
  }
  if (checkpoint.transformState.transformIdentifier !=
          layout.transformIdentifier() ||
      checkpoint.transformState.coefficientFamilies.size() !=
          layout.coefficientFamilyCount())
    return invalidShape(
        "Checkpoint template and output plan transform layouts differ.");
  for (std::size_t family = 0; family < layout.coefficientFamilyCount();
       ++family) {
    const auto &expected = layout.coefficientFamilies()[family];
    const auto &actual = checkpoint.transformState.coefficientFamilies[family];
    if (actual.identifier != expected.identifier ||
        actual.spectralDimensions != expected.spectralDimensions ||
        actual.values.size() != expected.elementCount)
      return invalidShape(
          "Checkpoint template coefficient families do not match the output "
          "plan layout.");
  }
  return WVKernelStatus::ok();
}

WVKernelStatus stageOutputCheckpointState(
    const WVIntegrationStateLayout &layout, WVCheckpoint &checkpoint,
    const WVIntegrationState &state) {
  auto status = validateOutputCheckpointTemplate(layout, checkpoint);
  if (!status)
    return status;
  if (!layout.hasLegacyCoefficientTriple()) {
    status = captureTransformStateCheckpoint(layout, state,
                                             checkpoint.transformState);
    if (!status)
      return status;
    checkpoint.state.t = state.waveVortex.t;
    checkpoint.state.t0 = state.waveVortex.t0;
    return WVKernelStatus::ok();
  }
  const auto count = layout.coefficientShape().elementCount();
  const WVComplexConstView sources[] = {state.waveVortex.coefficients.Ap,
                                        state.waveVortex.coefficients.Am,
                                        state.waveVortex.coefficients.A0};
  std::vector<WVComplex64> *destinations[] = {
      &checkpoint.state.coefficients.Ap, &checkpoint.state.coefficients.Am,
      &checkpoint.state.coefficients.A0};
  for (std::size_t component = 0; component < 3; ++component)
    std::copy_n(sources[component].data, count,
                destinations[component]->data());
  checkpoint.state.t = state.waveVortex.t;
  checkpoint.state.t0 = state.waveVortex.t0;
  return WVKernelStatus::ok();
}

} // namespace wavevortex::runtime::detail
