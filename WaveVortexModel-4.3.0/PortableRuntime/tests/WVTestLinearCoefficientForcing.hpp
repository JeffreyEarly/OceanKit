#pragma once

#include "WaveVortexRuntime/WVExtensionCatalog.hpp"

#include <cmath>
#include <memory>
#include <utility>

namespace wavevortex::runtime::test {

inline constexpr const char* LinearCoefficientForcingIdentifier =
    "WVTestPortableLinearCoefficientForcing";

class WVTestPortableLinearCoefficientForcing final : public WVForcing {
public:
  WVTestPortableLinearCoefficientForcing(WVFrozenForcingEntry entry,
                                         double rate)
      : typeIdentifier_(std::move(entry.typeIdentifier)),
        name_(std::move(entry.name)), contractVersion_(entry.contractVersion),
        stage_(entry.stage), priority_(entry.priority),
        ordinal_(entry.ordinal), rate_(rate) {}

  const std::string &typeIdentifier() const noexcept override {
    return typeIdentifier_;
  }
  std::uint32_t contractVersion() const noexcept override {
    return contractVersion_;
  }
  const std::string &name() const noexcept override { return name_; }
  WVForcingStage stage() const noexcept override { return stage_; }
  std::uint8_t priority() const noexcept override { return priority_; }
  std::size_t ordinal() const noexcept override { return ordinal_; }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + typeIdentifier_.capacity() + name_.capacity();
  }
  WVKernelStatus
  addRightHandSide(WVForcingExecutionContext &context) const override {
    return context.linearCoefficientTendency(rate_);
  }

private:
  std::string typeIdentifier_;
  std::string name_;
  std::uint32_t contractVersion_ = 0;
  WVForcingStage stage_ = WVForcingStage::spectral;
  std::uint8_t priority_ = 255;
  std::size_t ordinal_ = 0;
  double rate_ = 0.0;
};

inline WVKernelStatus createLinearCoefficientForcing(
    const WVFrozenForcingEntry &entry,
    const WVTransformConstantStratificationDescriptor &, bool,
    std::unique_ptr<WVForcing> &forcing) {
  const auto *value = entry.configuration.value("rate");
  const auto *values =
      value == nullptr
          ? nullptr
          : std::get_if<std::vector<double>>(&value->storage);
  if (values == nullptr || values->size() != 1 ||
      !std::isfinite(values->front()))
    return {WVKernelStatusCode::invalidConfiguration,
            "Test linear forcing requires one finite scalar rate."};
  forcing = std::make_unique<WVTestPortableLinearCoefficientForcing>(
      entry, values->front());
  return WVKernelStatus::ok();
}

inline WVForcingFactoryRegistration
linearCoefficientRegistration() {
  WVForcingPersistenceSchema persistence;
  persistence.writesNameAttribute = true;
  persistence.fields.push_back(
      {WVForcingPersistenceEncoding::realVariable, "rate", {}, "rate",
       WVForcingDimensionRule::scalar, {}, false, false});
  return {LinearCoefficientForcingIdentifier,
          WVPortablePairContractVersion,
          {"Spectral"},
          "test linear coefficient forcing",
          WVForcingStage::spectral,
          90,
          std::move(persistence),
          createLinearCoefficientForcing,
          true,
          "",
          false,
          WVForcingStage::spatial,
          {},
          {}};
}

inline WVKernelStatus registerLinearCoefficientForcing(
    WVExtensionCatalogBuilder &builder) {
  return builder.addForcingFactory(linearCoefficientRegistration());
}

} // namespace wavevortex::runtime::test
