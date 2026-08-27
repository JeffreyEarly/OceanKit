#pragma once

#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WVTestLinearCoefficientForcing.hpp"
#include "WVTestQuadraticSchedule.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>

namespace wavevortex::runtime::test {

inline std::shared_ptr<const WVExtensionCatalog> makeExtensionCatalog() {
  WVExtensionCatalogBuilder builder;
  auto status = addBuiltInExtensions(builder);
  if (!status)
    throw std::runtime_error("Unable to add built-in test extensions: " +
                             status.message);
  auto registrations = builtInForcingFactories();
  const auto fixed = std::find_if(
      registrations.begin(), registrations.end(), [](const auto &value) {
        return value.matlabClassName == "WVFixedAmplitudeForcing";
      });
  if (fixed == registrations.end())
    throw std::runtime_error("Built-in fixed-amplitude forcing is absent.");
  auto alias = *fixed;
  alias.matlabClassName = "WVTestPortableFixedAmplitudeForcing";
  status = builder.addForcingFactory(std::move(alias));
  if (status)
    status = registerLinearCoefficientForcing(builder);
  if (status)
    status = registerQuadraticSchedule(builder);
  if (!status)
    throw std::runtime_error("Unable to add test extensions: " +
                             status.message);
  std::shared_ptr<const WVExtensionCatalog> catalog;
  status = builder.freeze(catalog);
  if (!status)
    throw std::runtime_error("Unable to build the test extension catalog: " +
                             status.message);
  return catalog;
}

inline const std::shared_ptr<const WVExtensionCatalog> &extensionCatalog() {
  static const auto catalog = makeExtensionCatalog();
  return catalog;
}

} // namespace wavevortex::runtime::test
