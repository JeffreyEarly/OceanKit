#pragma once

#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace wavevortex::runtime::detail {

// Construction-resolved transform persistence boundary. The transactional
// sink consumes only ordered coefficient-family layouts and delegates root
// schema/configuration behavior to one named transform adapter.
class WVModelOutputTransformAdapter {
public:
  virtual ~WVModelOutputTransformAdapter() = default;
  virtual WVCheckpointStatus validate(
      const WVCheckpoint &checkpoint,
      const WVIntegrationStateLayout &layout) const = 0;
  virtual const WVComplex64 *coefficientData(
      const WVCheckpoint &checkpoint,
      const WVIntegrationStateLayout &layout,
      std::size_t family) const noexcept = 0;
  virtual void bindConstructionState(
      const WVCheckpoint &checkpoint,
      WVIntegrationState &state) const noexcept = 0;
  virtual WVCheckpointStatus defineRoot(
      int root, const WVCheckpoint &checkpoint,
      const WVForcingCatalog &catalog, bool isDynamicsLinear,
      std::array<int, 5> &dimensions, std::vector<int> &forcingGroups,
      std::vector<const WVFrozenForcingEntry *> &forcingEntries) const = 0;
  virtual WVCheckpointStatus writeRoot(
      int root, const WVCheckpoint &checkpoint,
      const WVForcingCatalog &catalog,
      const std::vector<int> &forcingGroups,
      const std::vector<const WVFrozenForcingEntry *> &forcingEntries) const = 0;
  virtual bool sameConfiguration(
      const WVCheckpointInspection &inspection) const noexcept = 0;
  virtual std::size_t persistentBytes() const noexcept = 0;
};

WVCheckpointStatus createModelOutputTransformAdapter(
    const WVCheckpoint &checkpoint,
    std::unique_ptr<WVModelOutputTransformAdapter> &adapter);
bool sameModelOutputTransformConfiguration(
    const WVCheckpointInspection &left,
    const WVCheckpointInspection &right) noexcept;
bool modelOutputGroupCarriesCompleteCoefficientRestart(
    const WVTransformStateDescription &description,
    bool hasDeclaredCoefficientFamilies,
    bool hasCoefficientObserver) noexcept;

WVCheckpointStatus defineModelOutputRoot(
    int root, const WVCheckpoint &checkpoint, const WVForcingCatalog &catalog,
    bool isDynamicsLinear,
    std::array<int, 5> &dimensions, std::vector<int> &forcingGroups,
    std::vector<const WVFrozenForcingEntry *> &forcingEntries);

WVCheckpointStatus writeModelOutputRoot(
    int root, const WVCheckpoint &checkpoint,
    const WVForcingCatalog &catalog,
    const std::vector<int> &forcingGroups,
    const std::vector<const WVFrozenForcingEntry *> &forcingEntries);

WVCheckpointStatus checkedNetCDF(int code, const std::string &operation,
                                 const std::string &location);
WVCheckpointStatus putTextAttribute(int group, int variable, const char *name,
                                    const std::string &value,
                                    const std::string &location);
WVCheckpointStatus putByteAttribute(int group, int variable, const char *name,
                                    unsigned char value,
                                    const std::string &location);
WVCheckpointStatus defineDoubleVariable(int group, const std::string &name,
                                        const std::vector<int> &dimensions,
                                        int &variable, const std::string &path);
WVCheckpointStatus defineComplexVariable(int group, const std::string &name,
                                         const std::vector<int> &dimensions,
                                         const std::string &path);

} // namespace wavevortex::runtime::detail
