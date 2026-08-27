#include "WaveVortexRuntime/WVModelOutputNetCDF.hpp"
#include "WaveVortexRuntime/generated/WVPortableVariableCatalog.hpp"

#include "WVLegacyObservationNetCDFAdapter.hpp"
#include "WVModelOutputTransformAdapters.hpp"
#include "WVNetCDF.hpp"

#include <netcdf.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <new>
#include <set>
#include <string>
#include <system_error>
#include <utility>

namespace wavevortex::runtime {
namespace {

WVCheckpointStatus failure(WVCheckpointStatusCode code, std::string message,
                           std::string location) {
  return {code, std::move(message), std::move(location)};
}

WVCheckpointStatus validateCatalogIdentity(
    const WVModelOutputNetCDFConfiguration &configuration,
    const WVPortableObserverDescriptor &descriptor) {
  if (!configuration.catalog)
    return failure(WVCheckpointStatusCode::schemaMismatch,
                   "Output configuration has no extension catalog.", "/");
  if (descriptor.catalog() != configuration.catalog)
    return failure(WVCheckpointStatusCode::schemaMismatch,
                   "The output configuration and observer descriptor require "
                   "the same extension catalog.",
                   "/");
  return WVCheckpointStatus::ok();
}

const WVPortableVariableMetadata *
coefficientMetadata(std::string_view name) noexcept {
  const auto *metadata = findPortableVariable(name);
  return metadata != nullptr &&
                 metadata->kind == WVPortableVariableKind::coefficient
             ? metadata
             : nullptr;
}

const char *coordinateRoleName(WVObservationCoordinateRole role) noexcept {
  switch (role) {
  case WVObservationCoordinateRole::none:
    return "none";
  case WVObservationCoordinateRole::recordTime:
    return "record-time";
  case WVObservationCoordinateRole::sampleTime:
    return "sample-time";
  case WVObservationCoordinateRole::x:
    return "x";
  case WVObservationCoordinateRole::y:
    return "y";
  case WVObservationCoordinateRole::z:
    return "z";
  case WVObservationCoordinateRole::identifier:
    return "identifier";
  case WVObservationCoordinateRole::depth:
    return "depth";
  case WVObservationCoordinateRole::pass:
    return "pass";
  case WVObservationCoordinateRole::profile:
    return "profile";
  }
  return "none";
}

const char *raggedRoleName(WVObservationRaggedRole role) noexcept {
  switch (role) {
  case WVObservationRaggedRole::none:
    return "none";
  case WVObservationRaggedRole::rowCount:
    return "row-count";
  case WVObservationRaggedRole::rowOffset:
    return "row-offset";
  }
  return "none";
}

const char *valueLayoutName(WVObservationValueLayout layout) noexcept {
  switch (layout) {
  case WVObservationValueLayout::staticValue:
    return "static";
  case WVObservationValueLayout::initialValue:
    return "initial";
  case WVObservationValueLayout::record:
    return "record";
  case WVObservationValueLayout::flat:
    return "flat";
  }
  return "record";
}

WVKernelStatus kernelFailure(const WVCheckpointStatus &status) {
  return {WVKernelStatusCode::unsupportedOperation,
          status.message + (status.location.empty()
                                ? std::string{}
                                : " at " + status.location)};
}

std::filesystem::path temporaryPath(const std::filesystem::path &destination) {
  static std::atomic<std::uint64_t> sequence{0};
  const auto stamp =
      std::chrono::steady_clock::now().time_since_epoch().count();
  return destination.parent_path() /
         (destination.filename().string() + ".wv-output-tmp-" +
          std::to_string(stamp) + "-" + std::to_string(++sequence));
}

WVCheckpointStatus
putStringListAttribute(int group, const char *name,
                       const std::vector<std::string> &values,
                       const std::string &path) {
  if (values.empty())
    return WVCheckpointStatus::ok();
  std::vector<const char *> raw;
  raw.reserve(values.size());
  for (const auto &value : values)
    raw.push_back(value.c_str());
  return detail::checkedNetCDF(
      nc_put_att_string(group, NC_GLOBAL, name, raw.size(), raw.data()),
      "String-list attribute definition", path + "/@" + name);
}

WVCheckpointStatus defineScalar(int group, const std::string &name,
                                int &variable, const std::string &path) {
  return detail::defineDoubleVariable(group, name, {}, variable, path);
}

WVCheckpointStatus writeScalar(int group, const std::string &name, double value,
                               const std::string &path) {
  int variable = -1;
  auto result =
      detail::checkedNetCDF(nc_inq_varid(group, name.c_str(), &variable),
                            "Variable lookup", path + "/" + name);
  if (!result)
    return result;
  return detail::checkedNetCDF(nc_put_var_double(group, variable, &value),
                               "Variable write", path + "/" + name);
}

std::string hexEncode(const std::vector<std::uint8_t> &bytes) {
  static constexpr char digits[] = "0123456789abcdef";
  std::string result(bytes.size() * 2, '0');
  for (std::size_t index = 0; index < bytes.size(); ++index) {
    result[2 * index] = digits[bytes[index] >> 4U];
    result[2 * index + 1] = digits[bytes[index] & 0x0fU];
  }
  return result;
}

bool hexDecode(const std::string &text, std::vector<std::uint8_t> &bytes) {
  if (text.size() % 2 != 0)
    return false;
  auto digit = [](char value) -> int {
    if (value >= '0' && value <= '9')
      return value - '0';
    if (value >= 'a' && value <= 'f')
      return value - 'a' + 10;
    if (value >= 'A' && value <= 'F')
      return value - 'A' + 10;
    return -1;
  };
  bytes.resize(text.size() / 2);
  for (std::size_t index = 0; index < bytes.size(); ++index) {
    const auto high = digit(text[2 * index]);
    const auto low = digit(text[2 * index + 1]);
    if (high < 0 || low < 0)
      return false;
    bytes[index] = static_cast<std::uint8_t>((high << 4) | low);
  }
  return true;
}

bool sameTypedRecord(const WVPortableTypedRecord &left,
                     const WVPortableTypedRecord &right) noexcept {
  if (left.schemaIdentifier != right.schemaIdentifier ||
      left.schemaVersion != right.schemaVersion ||
      left.values.size() != right.values.size())
    return false;
  for (std::size_t index = 0; index < left.values.size(); ++index) {
    const auto &a = left.values[index];
    const auto &b = right.values[index];
    if (a.name != b.name || a.dimensions != b.dimensions ||
        a.storage != b.storage)
      return false;
  }
  return true;
}

bool sameDestinationProgress(
    const std::vector<WVOutputDestinationProgress> &left,
    const std::vector<WVOutputDestinationProgress> &right) noexcept {
  if (left.size() != right.size())
    return false;
  for (std::size_t index = 0; index < left.size(); ++index) {
    const auto &a = left[index];
    const auto &b = right[index];
    if (a.fileIdentifier != b.fileIdentifier ||
        a.groupIdentifier != b.groupIdentifier ||
        a.recordCount != b.recordCount ||
        a.hasCommittedTime != b.hasCommittedTime ||
        (a.hasCommittedTime && a.lastCommittedTime != b.lastCommittedTime) ||
        a.committedScheduleCursor.committedOrdinal !=
            b.committedScheduleCursor.committedOrdinal ||
        !sameTypedRecord(a.committedScheduleCursor.values,
                         b.committedScheduleCursor.values) ||
        a.unlimitedAxes.size() != b.unlimitedAxes.size())
      return false;
    for (std::size_t axis = 0; axis < a.unlimitedAxes.size(); ++axis)
      if (a.unlimitedAxes[axis].axisIdentifier !=
              b.unlimitedAxes[axis].axisIdentifier ||
          a.unlimitedAxes[axis].committedCount !=
              b.unlimitedAxes[axis].committedCount ||
          a.unlimitedAxes[axis].physicalCount !=
              b.unlimitedAxes[axis].physicalCount)
        return false;
  }
  return true;
}

} // namespace

class WVModelOutputNetCDFSink::Impl {
public:
  struct Variable {
    std::string observerIdentifier;
    std::string schemaIdentifier;
    std::uint32_t schemaVersion = WVObservationSchemaContractVersion;
    WVObservationVariable specification;
    std::vector<WVObservationAxis> axes;
    // Resolved once before delivery. Entries align with axes; fixed axes use
    // WVNoResolvedObservationVariable. The write vectors use NetCDF order and
    // retain their bounded capacity between occurrences.
    std::vector<std::size_t> unlimitedAxisIndices;
    std::vector<std::size_t> writeStart;
    std::vector<std::size_t> writeCount;
    int realId = -1;
    int imagId = -1;
  };

  struct ObserverSchema {
    ObserverSchema(std::string identifier, WVObservationSchema value)
        : observerIdentifier(std::move(identifier)), schema(std::move(value)) {}
    std::string observerIdentifier;
    WVObservationSchema schema;
    std::size_t observerOrdinal = 0;
    std::vector<std::vector<std::size_t>> variableAxisIndices;
    std::vector<std::size_t> persistenceVariableIndices;
    std::vector<std::size_t> unlimitedAxisIndices;
    std::vector<std::size_t> raggedChildAxisIndices;
    // Reused validation workspace. These vectors are configuration-sized and
    // never grow with the number of output occurrences.
    std::vector<std::size_t> preparedValueIndices;
    std::vector<std::size_t> preparedUnlimitedExtents;
    std::vector<std::uint8_t> preparedObservedVariables;
  };

  struct UnlimitedAxis {
    std::string name;
    std::size_t committedCount = 0;
    std::size_t physicalCount = 0;
    int dimensionId = -1;
    int progressVariableId = -1;
  };

  struct Group {
    explicit Group(const WVOutputGroupRecord &canonicalRecord)
        : record(canonicalRecord) {}

    const WVOutputGroupRecord &record;
    int id = -1;
    int timeId = -1;
    int scheduleOrdinalId = -1;
    int scheduleCursorId = -1;
    std::vector<int> coefficientReal;
    std::vector<int> coefficientImag;
    std::vector<Variable> derivedVariables;
    std::vector<ObserverSchema> observerSchemas;
    std::vector<UnlimitedAxis> unlimitedAxes;
    std::size_t recordCount = 0;
    bool hasCommittedTime = false;
    double lastCommittedTime = 0.0;
    WVOutputScheduleOrdinal committedOrdinal = WVNoCommittedOutputOrdinal;
    WVPortableTypedRecord scheduleCursor;
    // Event-only numeric bindings compiled while all occurrence batches are
    // prepared, before this group can mutate its destination.
    std::vector<std::size_t> preparedBatchIndices;
    std::vector<std::size_t> preparedUnlimitedAxisIncrements;
    std::vector<double> realWriteScratch;
    std::vector<double> imaginaryWriteScratch;
    std::vector<long long> integerWriteScratch;
    std::vector<const char *> textWriteScratch;
  };

  struct File {
    explicit File(const WVOutputFileRecord &canonicalRecord)
        : record(canonicalRecord) {}

    const WVOutputFileRecord &record;
    std::filesystem::path destination;
    std::filesystem::path staging;
    int id = -1;
    std::vector<Group> groups;
  };

  struct PreparedObservationBatch {
    WVObservationOccurrenceIdentity identity;
    WVObservationBatch batch;
  };

  std::shared_ptr<const WVExtensionCatalog> catalog;
  std::unique_ptr<detail::WVModelOutputTransformAdapter> transformAdapter;
  const WVCheckpoint *constructionCheckpoint = nullptr;
  bool isDynamicsLinear = false;
  WVPortableObserverDescriptor descriptor;
  WVIntegrationStateLayout stateLayout;
  WVObserverSampleSource *sampleSource = nullptr;
  std::vector<File> files;
  std::vector<WVOutputDestinationProgress> destinationProgress;
  WVModelOutputNetCDFMetrics metrics;
  std::size_t preparedEventOrdinal = std::numeric_limits<std::size_t>::max();
  double preparedScheduledTime =
      std::numeric_limits<double>::quiet_NaN();
  std::size_t preparedEventRouteCount = 0;
  std::vector<std::pair<std::size_t, std::size_t>> committedPreparedRoutes;
  std::vector<PreparedObservationBatch> preparedObservationBatches;
  WVCheckpointStatus preparedObservationStatus = WVCheckpointStatus::ok();
  bool preparedEventComplete = false;
  bool appendMode = false;
  bool closed = false;
  bool preflightComplete = false;

  std::vector<std::string> coefficientFamilyIdentifiers() const {
    std::vector<std::string> identifiers;
    identifiers.reserve(stateLayout.coefficientFamilyCount());
    for (const auto &family : stateLayout.coefficientFamilies())
      identifiers.push_back(family.identifier);
    return identifiers;
  }

  const WVComplex64 *constructionCoefficientData(
      std::size_t family) const noexcept {
    if (constructionCheckpoint == nullptr ||
        transformAdapter == nullptr ||
        family >= stateLayout.coefficientFamilyCount())
      return nullptr;
    return transformAdapter->coefficientData(
        *constructionCheckpoint, stateLayout, family);
  }

  std::size_t sinkOccurrenceWorkspaceRetainedBytes() const noexcept {
    std::size_t bytes =
        preparedObservationBatches.capacity() *
        sizeof(PreparedObservationBatch);
    for (const auto &prepared : preparedObservationBatches)
      bytes += prepared.batch.metrics().retainedStorageBytes;
    for (const auto &file : files)
      for (const auto &group : file.groups)
        bytes += group.realWriteScratch.capacity() * sizeof(double) +
                 group.imaginaryWriteScratch.capacity() * sizeof(double) +
                 group.integerWriteScratch.capacity() * sizeof(long long) +
                 group.textWriteScratch.capacity() * sizeof(const char *);
    return bytes;
  }

  void updatePreparedObservationMetrics() {
    std::size_t retainedBytes = preparedObservationBatches.capacity() *
                                sizeof(PreparedObservationBatch);
    std::size_t liveBytes = 0;
    for (const auto &prepared : preparedObservationBatches) {
      const auto batchMetrics = prepared.batch.metrics();
      retainedBytes += batchMetrics.retainedStorageBytes;
      liveBytes += batchMetrics.liveBytes;
    }
    metrics.batchRetainedStorageBytes =
        std::max(metrics.batchRetainedStorageBytes, retainedBytes);
    metrics.batchMaximumLiveBytes =
        std::max(metrics.batchMaximumLiveBytes, liveBytes);
    const auto evaluatorRetained =
        sampleSource == nullptr
            ? 0
            : sampleSource->occurrenceWorkspaceRetainedBytes();
    const auto evaluatorLive =
        sampleSource == nullptr ? 0
                                : sampleSource->occurrenceWorkspaceLiveBytes();
    const auto sinkWorkspaceRetained =
        sinkOccurrenceWorkspaceRetainedBytes();
    metrics.occurrenceWorkspaceRetainedBytes =
        evaluatorRetained + sinkWorkspaceRetained;
    metrics.occurrenceWorkspaceMaximumLiveBytes = std::max(
        metrics.occurrenceWorkspaceMaximumLiveBytes,
        std::max(evaluatorRetained, evaluatorLive) +
            sinkWorkspaceRetained);
    updateRetainedStorageMetric();
  }

  static const WVObservationAxis *
  axis(const WVObservationSchema &schema,
       const std::string &identifier) noexcept {
    const auto found =
        std::find_if(schema.axes.begin(), schema.axes.end(),
                     [&](const auto &candidate) {
                       return candidate.identifier == identifier;
                     });
    return found == schema.axes.end() ? nullptr : &*found;
  }

  static bool applicable(const WVObservationVariable &variable,
                         WVObservationBatchKind kind) noexcept {
    const bool initial =
        variable.layout == WVObservationValueLayout::staticValue ||
        variable.layout == WVObservationValueLayout::initialValue;
    return kind == WVObservationBatchKind::initial ? initial : !initial;
  }

  WVKernelStatus compileVariableWriteBindings(Group &group) {
    constexpr auto unresolved = WVNoResolvedObservationVariable;
    for (auto &variable : group.derivedVariables) {
      variable.unlimitedAxisIndices.assign(variable.axes.size(), unresolved);
      for (std::size_t axisIndex = 0; axisIndex < variable.axes.size();
           ++axisIndex) {
        const auto &axisDefinition = variable.axes[axisIndex];
        if (axisDefinition.kind != WVObservationAxisKind::unlimited)
          continue;
        const auto found = std::find_if(
            group.unlimitedAxes.begin(), group.unlimitedAxes.end(),
            [&](const auto &candidate) {
              return candidate.name == axisDefinition.name;
            });
        if (found == group.unlimitedAxes.end())
          return {WVKernelStatusCode::invalidConfiguration,
                  "An observation variable has no resolved unlimited axis."};
        variable.unlimitedAxisIndices[axisIndex] =
            static_cast<std::size_t>(found - group.unlimitedAxes.begin());
      }
      const auto persistedRank =
          variable.axes.size() +
          (variable.specification.layout == WVObservationValueLayout::record
               ? 1
               : 0);
      variable.writeStart.assign(persistedRank, 0);
      variable.writeCount.assign(persistedRank, 1);
    }
    return WVKernelStatus::ok();
  }

  WVKernelStatus compileEventBindings(Group &group,
                                      const WVOutputRouteView &route) {
    constexpr auto unresolved = WVNoResolvedObservationVariable;
    if (route.observerCount != group.observerSchemas.size())
      return {WVKernelStatusCode::invalidConfiguration,
              "Output route and observation schemas differ."};
    const auto writeBindingStatus = compileVariableWriteBindings(group);
    if (!writeBindingStatus)
      return writeBindingStatus;

    for (std::size_t observerIndex = 0; observerIndex < route.observerCount;
         ++observerIndex) {
      const auto &observerView = route.observers[observerIndex];
      if (observerView.record == nullptr)
        return {WVKernelStatusCode::invalidConfiguration,
                "Output route has no resolved observer."};
      auto &binding = group.observerSchemas[observerIndex];
      binding.observerOrdinal = observerView.observerOrdinal;
      binding.unlimitedAxisIndices.assign(binding.schema.axes.size(),
                                           unresolved);
      for (std::size_t axisIndex = 0; axisIndex < binding.schema.axes.size();
           ++axisIndex) {
        const auto &axisDefinition = binding.schema.axes[axisIndex];
        if (axisDefinition.kind != WVObservationAxisKind::unlimited)
          continue;
        const auto found = std::find_if(
            group.unlimitedAxes.begin(), group.unlimitedAxes.end(),
            [&](const auto &candidate) {
              return candidate.name == axisDefinition.name;
            });
        if (found == group.unlimitedAxes.end())
          return {WVKernelStatusCode::invalidConfiguration,
                  "Observer schema has no resolved unlimited axis."};
        binding.unlimitedAxisIndices[axisIndex] =
            static_cast<std::size_t>(found - group.unlimitedAxes.begin());
      }

      binding.variableAxisIndices.assign(binding.schema.variables.size(), {});
      binding.persistenceVariableIndices.assign(binding.schema.variables.size(),
                                                unresolved);
      binding.raggedChildAxisIndices.assign(binding.schema.variables.size(),
                                            unresolved);
      for (std::size_t variableIndex = 0;
           variableIndex < binding.schema.variables.size(); ++variableIndex) {
        const auto &specification = binding.schema.variables[variableIndex];
        auto &axisIndices = binding.variableAxisIndices[variableIndex];
        axisIndices.reserve(specification.dimensionIdentifiers.size());
        for (const auto &identifier : specification.dimensionIdentifiers) {
          const auto found = std::find_if(
              binding.schema.axes.begin(), binding.schema.axes.end(),
              [&](const auto &candidate) {
                return candidate.identifier == identifier;
              });
          if (found == binding.schema.axes.end())
            return {WVKernelStatusCode::invalidConfiguration,
                    "Observer variable has an unresolved axis."};
          axisIndices.push_back(
              static_cast<std::size_t>(found - binding.schema.axes.begin()));
        }
        if (specification.raggedRole != WVObservationRaggedRole::none) {
          const auto found = std::find_if(
              binding.schema.axes.begin(), binding.schema.axes.end(),
              [&](const auto &candidate) {
                return candidate.identifier ==
                       specification.raggedChildAxisIdentifier;
              });
          if (found == binding.schema.axes.end())
            return {WVKernelStatusCode::invalidConfiguration,
                    "Ragged observation variable has an unresolved child axis."};
          binding.raggedChildAxisIndices[variableIndex] =
              static_cast<std::size_t>(found - binding.schema.axes.begin());
        }
        const auto persisted = std::find_if(
            group.derivedVariables.begin(), group.derivedVariables.end(),
            [&](const auto &candidate) {
              return candidate.observerIdentifier ==
                         binding.observerIdentifier &&
                     candidate.specification.identifier ==
                         specification.identifier;
            });
        if (persisted != group.derivedVariables.end())
          binding.persistenceVariableIndices[variableIndex] =
              static_cast<std::size_t>(persisted -
                                       group.derivedVariables.begin());
        else if (!(group.record.containsCompleteCoefficientRestart &&
                   coefficientMetadata(specification.name) != nullptr))
          return {WVKernelStatusCode::invalidConfiguration,
                  "Observer variable has no resolved persistence binding."};
      }
      binding.preparedValueIndices.assign(binding.schema.variables.size(),
                                          unresolved);
      binding.preparedUnlimitedExtents.assign(binding.schema.axes.size(),
                                              unresolved);
      binding.preparedObservedVariables.assign(binding.schema.variables.size(),
                                               0);
    }
    group.preparedUnlimitedAxisIncrements.assign(group.unlimitedAxes.size(),
                                                 unresolved);
    return WVKernelStatus::ok();
  }

  WVKernelStatus prepareEventBindingContracts(Group &group) {
    if (sampleSource == nullptr)
      return WVKernelStatus::ok();
    for (const auto &identifier : group.record.observerIdentifiers) {
      const auto *record = observer(identifier);
      if (record == nullptr)
        return {WVKernelStatusCode::invalidConfiguration,
                "An output group references an unknown observer."};
      WVObservationSchema schema;
      auto status = sampleSource->observationSchema(*record, schema);
      if (!status)
        return status;
      status = validateObservationSchema(schema);
      if (!status)
        return status;

      for (const auto &axisDefinition : schema.axes) {
        if (axisDefinition.kind != WVObservationAxisKind::unlimited)
          continue;
        const auto found =
            std::find_if(group.unlimitedAxes.begin(), group.unlimitedAxes.end(),
                         [&](const auto &candidate) {
                           return candidate.name == axisDefinition.name;
                         });
        if (found == group.unlimitedAxes.end())
          group.unlimitedAxes.push_back({axisDefinition.name});
      }
      for (const auto &specification : schema.variables) {
        if (group.record.containsCompleteCoefficientRestart &&
            coefficientMetadata(specification.name) != nullptr)
          continue;
        Variable variable;
        variable.observerIdentifier = identifier;
        variable.schemaIdentifier = schema.identifier;
        variable.schemaVersion = schema.version;
        variable.specification = specification;
        variable.axes.reserve(specification.dimensionIdentifiers.size());
        for (const auto &axisIdentifier : specification.dimensionIdentifiers) {
          const auto *axisDefinition = axis(schema, axisIdentifier);
          if (axisDefinition == nullptr)
            return {WVKernelStatusCode::invalidConfiguration,
                    "An observation variable has an unresolved axis."};
          variable.axes.push_back(*axisDefinition);
        }
        group.derivedVariables.push_back(std::move(variable));
      }
      group.observerSchemas.emplace_back(identifier, std::move(schema));
    }
    return WVKernelStatus::ok();
  }

  WVKernelStatus validateAndCompilePlan(const WVOutputPlan &plan) {
    const auto &outputFiles = descriptor.record().outputFiles;
    if (plan.metrics().fileCount != outputFiles.size())
      return {WVKernelStatusCode::invalidConfiguration,
              "Output plan and NetCDF file graph differ."};
    std::size_t expectedGroupCount = 0;
    for (const auto &file : outputFiles)
      expectedGroupCount += file.groups.size();
    if (plan.groupCount() != expectedGroupCount)
      return {WVKernelStatusCode::invalidConfiguration,
              "Output plan and NetCDF group graph differ."};
    if (!files.empty()) {
      if (files.size() != outputFiles.size())
        return {WVKernelStatusCode::invalidConfiguration,
                "Output plan and prepared NetCDF file graph differ."};
      for (std::size_t fileIndex = 0; fileIndex < files.size(); ++fileIndex)
        if (files[fileIndex].groups.size() !=
            outputFiles[fileIndex].groups.size())
          return {WVKernelStatusCode::invalidConfiguration,
                  "Output plan and prepared NetCDF group graph differ."};
    }

    std::vector<std::vector<std::uint8_t>> visited;
    visited.reserve(outputFiles.size());
    for (const auto &file : outputFiles)
      visited.emplace_back(file.groups.size(), 0);
    for (std::size_t groupIndex = 0; groupIndex < plan.groupCount();
         ++groupIndex) {
      const auto route = plan.groupRoute(groupIndex);
      if (route.fileOrdinal >= outputFiles.size() ||
          route.groupOrdinal >= outputFiles[route.fileOrdinal].groups.size())
        return {WVKernelStatusCode::invalidConfiguration,
                "Output plan contains an out-of-range NetCDF route."};
      if (visited[route.fileOrdinal][route.groupOrdinal] != 0)
        return {WVKernelStatusCode::invalidConfiguration,
                "Output plan contains a duplicate NetCDF route."};
      visited[route.fileOrdinal][route.groupOrdinal] = 1;
      const auto &fileRecord = outputFiles[route.fileOrdinal];
      const auto &groupRecord = fileRecord.groups[route.groupOrdinal];
      if (route.fileIdentifier != fileRecord.identifier ||
          route.destination != fileRecord.destination ||
          route.groupIdentifier != groupRecord.identifier ||
          route.groupName != groupRecord.name ||
          route.observerCount != groupRecord.observerIdentifiers.size())
        return {WVKernelStatusCode::invalidConfiguration,
                "Output plan and NetCDF route identity differ."};
      for (std::size_t observerIndex = 0; observerIndex < route.observerCount;
           ++observerIndex) {
        if (route.observers[observerIndex].record == nullptr ||
            route.observers[observerIndex].record->identifier !=
                groupRecord.observerIdentifiers[observerIndex])
          return {WVKernelStatusCode::invalidConfiguration,
                  "Output plan and NetCDF observer route differ."};
      }
      if (sampleSource == nullptr)
        continue;
      if (files.empty()) {
        Group transient(groupRecord);
        const auto preparation = prepareEventBindingContracts(transient);
        if (!preparation)
          return preparation;
        const auto compilation = compileEventBindings(transient, route);
        if (!compilation)
          return compilation;
      } else {
        const auto compilation = compileEventBindings(
            files[route.fileOrdinal].groups[route.groupOrdinal], route);
        if (!compilation)
          return compilation;
      }
    }
    return WVKernelStatus::ok();
  }

  WVCheckpointStatus
  validateResolvedObservationBatch(ObserverSchema &binding,
                                   const WVObservationBatch &batch) {
    constexpr auto unresolved = WVNoResolvedObservationVariable;
    // The observer/schema association and every variable slot were compiled
    // before destination mutation. Event batches intentionally carry no
    // required names; numeric slots are the authoritative hot-path identity.
    if (batch.kind != WVObservationBatchKind::event ||
        batch.schemaVersion != binding.schema.version)
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     "Observation batch phase or schema version drifted.",
                     "/");
    std::size_t expectedValueCount = 0;
    for (const auto &variable : binding.schema.variables)
      expectedValueCount += applicable(variable, batch.kind) ? 1 : 0;
    if (batch.values.size() != expectedValueCount)
      return failure(
          WVCheckpointStatusCode::shapeMismatch,
          "Observation batch does not contain every applicable variable.",
          "/");
    std::fill(binding.preparedValueIndices.begin(),
              binding.preparedValueIndices.end(), unresolved);
    std::fill(binding.preparedUnlimitedExtents.begin(),
              binding.preparedUnlimitedExtents.end(), unresolved);
    std::fill(binding.preparedObservedVariables.begin(),
              binding.preparedObservedVariables.end(), 0);

    for (std::size_t valueIndex = 0; valueIndex < batch.values.size();
         ++valueIndex) {
      const auto &value = batch.values[valueIndex];
      const auto variableIndex = value.resolvedVariableIndex;
      if (variableIndex >= binding.schema.variables.size())
        return failure(WVCheckpointStatusCode::schemaMismatch,
                       "Observation value has no resolved variable slot.",
                       "/");
      const auto &variable = binding.schema.variables[variableIndex];
      if (binding.preparedObservedVariables[variableIndex] != 0 ||
          !applicable(variable, batch.kind))
        return failure(
            WVCheckpointStatusCode::schemaMismatch,
            "Observation value is duplicated or belongs to another phase.",
            "/");
      if (value.scalarType != variable.scalarType ||
          value.extents.size() !=
              binding.variableAxisIndices[variableIndex].size())
        return failure(WVCheckpointStatusCode::shapeMismatch,
                       "Observation value type or rank differs from its schema.",
                       "/");
      std::size_t expectedElements = 1;
      for (std::size_t dimension = 0; dimension < value.extents.size();
           ++dimension) {
        const auto axisIndex =
            binding.variableAxisIndices[variableIndex][dimension];
        const auto &axisDefinition = binding.schema.axes[axisIndex];
        const auto extent = value.extents[dimension];
        if (axisDefinition.kind == WVObservationAxisKind::fixed) {
          if (extent != axisDefinition.extent)
            return failure(WVCheckpointStatusCode::shapeMismatch,
                           "Observation fixed-axis extent changed.", "/");
        } else if (binding.preparedUnlimitedExtents[axisIndex] == unresolved) {
          binding.preparedUnlimitedExtents[axisIndex] = extent;
        } else if (binding.preparedUnlimitedExtents[axisIndex] != extent) {
          return failure(WVCheckpointStatusCode::shapeMismatch,
                         "Observation unlimited-axis extents disagree.", "/");
        }
        if (extent != 0 &&
            expectedElements >
                std::numeric_limits<std::size_t>::max() / extent)
          return failure(WVCheckpointStatusCode::invalidValue,
                         "Observation element count overflows size_t.", "/");
        expectedElements *= extent;
      }
      if (value.elementCount() != expectedElements)
        return failure(WVCheckpointStatusCode::invalidValue,
                       "Observation element count overflowed.", "/");
      const bool hasData =
          expectedElements == 0 ||
          (value.scalarType == WVObservationScalarType::real64 &&
           value.real64Data() != nullptr) ||
          (value.scalarType == WVObservationScalarType::complex64 &&
           value.complex64Data() != nullptr) ||
          (value.scalarType == WVObservationScalarType::integer64 &&
           value.integer64Data() != nullptr) ||
          (value.scalarType == WVObservationScalarType::boolean8 &&
           value.boolean8Data() != nullptr) ||
          (value.scalarType == WVObservationScalarType::text &&
           value.textData() != nullptr);
      if (!hasData)
        return failure(WVCheckpointStatusCode::invalidValue,
                       "Observation value has no compatible buffer.", "/");
      if (value.ownership == WVObservationBufferOwnership::owned) {
        const std::size_t ownedCount =
            value.scalarType == WVObservationScalarType::real64
                ? value.ownedReal64.size()
                : value.scalarType == WVObservationScalarType::complex64
                      ? value.ownedComplex64.size()
                      : value.scalarType == WVObservationScalarType::integer64
                            ? value.ownedInteger64.size()
                            : value.scalarType ==
                                      WVObservationScalarType::boolean8
                                  ? value.ownedBoolean8.size()
                                  : value.ownedText.size();
        if (ownedCount != expectedElements)
          return failure(WVCheckpointStatusCode::shapeMismatch,
                         "Owned observation buffer size differs from its extents.",
                         "/");
      }
      if (value.scalarType == WVObservationScalarType::boolean8) {
        const auto *booleans = value.boolean8Data();
        for (std::size_t index = 0; index < expectedElements; ++index)
          if (booleans[index] > 1)
            return failure(WVCheckpointStatusCode::invalidValue,
                           "Boolean observations must be zero or one.", "/");
      }
      binding.preparedObservedVariables[variableIndex] = 1;
      binding.preparedValueIndices[variableIndex] = valueIndex;
    }

    for (std::size_t variableIndex = 0;
         variableIndex < binding.schema.variables.size(); ++variableIndex) {
      const auto &variable = binding.schema.variables[variableIndex];
      if (!applicable(variable, batch.kind) ||
          variable.raggedRole == WVObservationRaggedRole::none)
        continue;
      const auto valueIndex = binding.preparedValueIndices[variableIndex];
      const auto childAxisIndex =
          binding.raggedChildAxisIndices[variableIndex];
      if (valueIndex >= batch.values.size() ||
          childAxisIndex >= binding.preparedUnlimitedExtents.size() ||
          binding.preparedUnlimitedExtents[childAxisIndex] == unresolved)
        return failure(WVCheckpointStatusCode::shapeMismatch,
                       "Ragged observation child axis is absent.", "/");
      const auto &value = batch.values[valueIndex];
      const auto *integers = value.integer64Data();
      const auto count = value.elementCount();
      const auto childExtent =
          binding.preparedUnlimitedExtents[childAxisIndex];
      if (variable.raggedRole == WVObservationRaggedRole::rowCount) {
        std::size_t total = 0;
        for (std::size_t index = 0; index < count; ++index) {
          if (integers[index] < 0 ||
              static_cast<std::uint64_t>(integers[index]) >
                  std::numeric_limits<std::size_t>::max())
            return failure(WVCheckpointStatusCode::invalidValue,
                           "Ragged row counts must be nonnegative and bounded.",
                           "/");
          const auto row = static_cast<std::size_t>(integers[index]);
          if (total > std::numeric_limits<std::size_t>::max() - row)
            return failure(WVCheckpointStatusCode::invalidValue,
                           "Ragged row-count sum overflows size_t.", "/");
          total += row;
        }
        if (total != childExtent)
          return failure(WVCheckpointStatusCode::shapeMismatch,
                         "Ragged row counts do not span their child axis.",
                         "/");
      } else {
        if (count == 0 && childExtent != 0)
          return failure(WVCheckpointStatusCode::shapeMismatch,
                         "An empty ragged parent has a nonempty child axis.",
                         "/");
        if (count != 0 && integers[0] != 0)
          return failure(WVCheckpointStatusCode::invalidValue,
                         "Ragged row offsets must begin at zero.", "/");
        for (std::size_t index = 0; index < count; ++index)
          if (integers[index] < 0 ||
              static_cast<std::uint64_t>(integers[index]) > childExtent ||
              (index > 0 && integers[index] < integers[index - 1]))
            return failure(WVCheckpointStatusCode::invalidValue,
                           "Ragged row offsets are malformed.", "/");
      }
    }
    return WVCheckpointStatus::ok();
  }

  WVCheckpointStatus prepareObservationBatches(const WVOutputEvent &event) {
    if (sampleSource == nullptr)
      return WVCheckpointStatus::ok();
    if (preparedEventOrdinal == event.eventOrdinal &&
        preparedScheduledTime == event.scheduledTime &&
        preparedEventRouteCount == event.routeCount) {
      if (preparedEventComplete)
        return failure(WVCheckpointStatusCode::appendConflict,
                       "A completed output event cannot be delivered twice.",
                       "/");
      return preparedObservationStatus;
    }

    preparedEventOrdinal = event.eventOrdinal;
    preparedScheduledTime = event.scheduledTime;
    preparedEventRouteCount = event.routeCount;
    preparedEventComplete = false;
    committedPreparedRoutes.clear();
    committedPreparedRoutes.reserve(event.routeCount);
    preparedObservationBatches.clear();
    preparedObservationStatus = WVCheckpointStatus::ok();
    const auto sourceStatus = sampleSource->prepare(event);
    if (!sourceStatus) {
      updatePreparedObservationMetrics();
      preparedObservationStatus = failure(
          WVCheckpointStatusCode::unsupportedObserver, sourceStatus.message,
          "/");
      return preparedObservationStatus;
    }

    std::size_t requestedBatchCount = 0;
    for (std::size_t routeIndex = 0; routeIndex < event.routeCount;
         ++routeIndex)
      requestedBatchCount += event.routes[routeIndex].observerCount;
    preparedObservationBatches.reserve(requestedBatchCount);
    for (std::size_t routeIndex = 0; routeIndex < event.routeCount;
         ++routeIndex) {
      const auto &route = event.routes[routeIndex];
      if (route.fileOrdinal >= files.size() ||
          route.groupOrdinal >= files[route.fileOrdinal].groups.size()) {
        preparedObservationStatus = failure(
            WVCheckpointStatusCode::schemaMismatch,
            "An event route is outside the configured output graph.", "/");
        return preparedObservationStatus;
      }
      auto &group = files[route.fileOrdinal].groups[route.groupOrdinal];
      if (route.observerCount != group.observerSchemas.size()) {
        preparedObservationStatus = failure(
            WVCheckpointStatusCode::schemaMismatch,
            "An event route has an incompatible observer schema count.", "/");
        return preparedObservationStatus;
      }
      group.preparedBatchIndices.assign(route.observerCount,
                                        std::numeric_limits<std::size_t>::max());
      for (std::size_t observerIndex = 0;
           observerIndex < route.observerCount; ++observerIndex) {
        const auto &observerView = route.observers[observerIndex];
        if (observerView.record == nullptr) {
          preparedObservationStatus = failure(
              WVCheckpointStatusCode::schemaMismatch,
              "An event route has no resolved observer record.", "/");
          return preparedObservationStatus;
        }
        WVObservationOccurrenceIdentity identity;
        const auto identityStatus = sampleSource->preparedOccurrenceIdentity(
            route, observerView, identity);
        if (!identityStatus) {
          updatePreparedObservationMetrics();
          preparedObservationStatus = failure(
              WVCheckpointStatusCode::unsupportedObserver,
              identityStatus.message, "/");
          return preparedObservationStatus;
        }
        if (identity.preparationOwner == nullptr ||
            identity.preparationGeneration == 0) {
          updatePreparedObservationMetrics();
          preparedObservationStatus = failure(
              WVCheckpointStatusCode::unsupportedObserver,
              "An observation source returned no exact prepared-occurrence "
              "cache token.", "/");
          return preparedObservationStatus;
        }
        std::size_t batchIndex = preparedObservationBatches.size();
        for (std::size_t candidate = 0;
             candidate < preparedObservationBatches.size(); ++candidate)
          if (samePreparedObservationOccurrenceIdentity(
                  preparedObservationBatches[candidate].identity, identity)) {
            batchIndex = candidate;
            break;
          }
        if (batchIndex == preparedObservationBatches.size()) {
          WVObservationBatch batch;
          const auto batchSourceStatus = sampleSource->observationBatch(
              identity, *observerView.record, batch);
          if (!batchSourceStatus) {
            updatePreparedObservationMetrics();
            preparedObservationStatus = failure(
                WVCheckpointStatusCode::unsupportedObserver,
                batchSourceStatus.message, "/");
            return preparedObservationStatus;
          }
          preparedObservationBatches.push_back(
              {identity, std::move(batch)});
        }
        const auto batchStatus = validateResolvedObservationBatch(
            group.observerSchemas[observerIndex],
            preparedObservationBatches[batchIndex].batch);
        if (!batchStatus) {
          updatePreparedObservationMetrics();
          preparedObservationStatus = batchStatus;
          return preparedObservationStatus;
        }
        group.preparedBatchIndices[observerIndex] = batchIndex;
      }
    }
    updatePreparedObservationMetrics();
    return preparedObservationStatus;
  }

  void completePreparedRoute(const WVOutputEvent &event,
                             const WVOutputRouteView &route) {
    if (sampleSource == nullptr)
      return;
    const auto routeKey =
        std::make_pair(route.fileOrdinal, route.groupOrdinal);
    if (std::find(committedPreparedRoutes.begin(),
                  committedPreparedRoutes.end(), routeKey) ==
        committedPreparedRoutes.end())
      committedPreparedRoutes.push_back(routeKey);
    if (committedPreparedRoutes.size() != preparedEventRouteCount)
      return;
    preparedObservationBatches.clear();
    committedPreparedRoutes.clear();
    preparedEventComplete = true;
    sampleSource->complete(event);
    updatePreparedObservationMetrics();
  }

  ~Impl() { close(); }

  const WVObserverRecord *observer(const std::string &identifier) const {
    const auto &record = descriptor.record();
    const auto found = std::find_if(record.observers.begin(),
                                    record.observers.end(),
                                    [&](const auto &candidate) {
                                      return candidate.identifier == identifier;
                                    });
    return found == record.observers.end() ? nullptr : &*found;
  }

  const WVStateBlockRecord *
  stateBlockRecord(const std::string &identifier) const {
    const auto &blocks = stateLayout.stateBlockRecords();
    const auto found =
        std::find_if(blocks.begin(), blocks.end(), [&](const auto &candidate) {
          return candidate.identifier == identifier;
        });
    return found == blocks.end() ? nullptr : &*found;
  }

  WVCheckpointStatus validateObservationGraph() const {
    if (sampleSource == nullptr)
      return WVCheckpointStatus::ok();
    std::map<std::string, WVObservationSchema> schemas;
    for (const auto &record : descriptor.record().observers) {
      WVObservationSchema schema;
      const auto sourceStatus = sampleSource->observationSchema(record, schema);
      if (!sourceStatus)
        return failure(WVCheckpointStatusCode::unsupportedObserver,
                       sourceStatus.message,
                       "/observingSystems/" + record.identifier);
      const auto schemaStatus = validateObservationSchema(schema);
      if (!schemaStatus)
        return failure(WVCheckpointStatusCode::schemaMismatch,
                       schemaStatus.message,
                       "/observingSystems/" + record.identifier);
      schemas.emplace(record.identifier, std::move(schema));
    }
    struct AxisContract {
      WVObservationAxisKind kind = WVObservationAxisKind::fixed;
      std::size_t extent = 0;
      WVObservationCoordinateRole role =
          WVObservationCoordinateRole::none;
    };
    for (const auto &file : descriptor.record().outputFiles)
        for (const auto &group : file.groups) {
        std::map<std::string, AxisContract> axesByPersistedName;
        std::map<std::string, std::string> variableOwners;
        variableOwners.emplace("t", "runtime");
        for (const auto &observerIdentifier : group.observerIdentifiers) {
          const auto schema = schemas.find(observerIdentifier);
          if (schema == schemas.end())
            return failure(WVCheckpointStatusCode::schemaMismatch,
                           "Output group references an unknown observer schema.",
                           file.destination + "/" + group.name);
          for (const auto &axisDefinition : schema->second.axes) {
            const AxisContract contract{axisDefinition.kind,
                                        axisDefinition.extent,
                                        axisDefinition.coordinateRole};
            const auto inserted = axesByPersistedName.emplace(
                axisDefinition.name, contract);
            if (!inserted.second &&
                (inserted.first->second.kind != contract.kind ||
                 inserted.first->second.extent != contract.extent))
              return failure(
                  WVCheckpointStatusCode::schemaMismatch,
                  "Observers declare incompatible axes with the same "
                  "persisted name.",
                  file.destination + "/" + group.name + "/" +
                      axisDefinition.name);
          }
          for (const auto &variable : schema->second.variables) {
            const bool canonicalCoefficient =
                group.containsCompleteCoefficientRestart &&
                coefficientMetadata(variable.name) != nullptr;
            if (canonicalCoefficient)
              continue;
            const auto claim = [&](const std::string &persistedName) {
              const auto inserted = variableOwners.emplace(
                  persistedName, observerIdentifier);
              if (inserted.second)
                return WVCheckpointStatus::ok();
              return failure(
                  WVCheckpointStatusCode::schemaMismatch,
                  "Observers declare duplicate persisted variable names "
                  "without an alias contract.",
                  file.destination + "/" + group.name + "/" +
                      persistedName);
            };
            if (variable.scalarType ==
                WVObservationScalarType::complex64) {
              auto result = claim(variable.name + "_real");
              if (!result)
                return result;
              result = claim(variable.name + "_imag");
              if (!result)
                return result;
            } else {
              const auto result = claim(variable.name);
              if (!result)
                return result;
            }
          }
        }
      }
    return WVCheckpointStatus::ok();
  }

  WVCheckpointStatus validateConfiguration() const {
    if (!catalog)
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     "Output configuration has no extension catalog.", "/");
    if (constructionCheckpoint == nullptr)
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     "Output construction has no checkpoint template.", "/");
    if (constructionCheckpoint->metadata.profileIdentifier !=
            WVCheckpointProfileIdentifier ||
        constructionCheckpoint->metadata.profileVersion !=
            WVCheckpointProfileVersion)
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     "Unsupported checkpoint template profile.", "/");
    if (transformAdapter == nullptr)
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     "Output construction has no resolved transform adapter.",
                     "/");
    auto transformStatus =
        transformAdapter->validate(*constructionCheckpoint, stateLayout);
    if (!transformStatus)
      return transformStatus;
    std::vector<double> coefficientTolerances;
    coefficientTolerances.reserve(stateLayout.coefficientFamilyCount());
    std::size_t coefficientIndex = 0;
    for (const auto &block : stateLayout.stateBlockRecords()) {
      const auto family = std::find_if(
          stateLayout.coefficientFamilies().begin(),
          stateLayout.coefficientFamilies().end(), [&](const auto &candidate) {
            return candidate.identifier == block.identifier;
          });
      if (family != stateLayout.coefficientFamilies().end()) {
        if (coefficientIndex >= stateLayout.coefficientFamilyCount() ||
            block.toleranceKind != WVToleranceKind::coefficientEnergyScaled ||
            !std::isfinite(block.absoluteTolerance) ||
            block.absoluteTolerance <= 0.0)
          return failure(
              WVCheckpointStatusCode::schemaMismatch,
              "Coefficient state blocks must carry a finite positive "
              "energy-scaled absTolerance for MATLAB persistence.",
              "/observingSystems/WVCoefficients/absTolerance");
        coefficientTolerances.push_back(block.absoluteTolerance);
        ++coefficientIndex;
      }
    }
    if (!isDynamicsLinear &&
        (coefficientIndex != stateLayout.coefficientFamilyCount() ||
         !std::all_of(coefficientTolerances.begin(),
                      coefficientTolerances.end(), [&](double value) {
                        return value == coefficientTolerances.front();
                      })))
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     "Coefficient families must use one common absTolerance.",
                     "/observingSystems/WVCoefficients/absTolerance");
    const auto coefficientFamilies = coefficientFamilyIdentifiers();
    for (const auto &record : descriptor.record().observers) {
      const bool canonicalCoefficientObserver =
          record.stateBlockIdentifiers == coefficientFamilies &&
          record.fieldNames.empty() && record.x.empty() && record.y.empty() &&
          record.z.empty();
      if (!canonicalCoefficientObserver && sampleSource == nullptr)
        return failure(WVCheckpointStatusCode::unsupportedObserver,
                       "Observer " + record.identifier +
                           " requires an observation sample source.",
                       "/observingSystems");
    }
    for (const auto &file : descriptor.record().outputFiles) {
      const auto restart = std::find_if(
          file.groups.begin(), file.groups.end(), [](const auto &group) {
            return group.containsCompleteCoefficientRestart;
          });
      if (restart == file.groups.end())
        return failure(WVCheckpointStatusCode::schemaMismatch,
                       "Each output file requires one complete restart group.",
                       file.destination);
      std::set<std::string> restartBlocks(coefficientFamilies.begin(),
                                          coefficientFamilies.end());
      for (const auto &group : file.groups) {
        for (const auto &observerIdentifier : group.observerIdentifiers) {
          const auto *record = observer(observerIdentifier);
          if (record != nullptr)
            restartBlocks.insert(record->stateBlockIdentifiers.begin(),
                                 record->stateBlockIdentifiers.end());
        }
      }
      for (const auto &block : stateLayout.stateBlockRecords()) {
        if (block.restartRequirement ==
                WVRestartRequirement::requiredDynamicState &&
            restartBlocks.find(block.identifier) == restartBlocks.end())
          return failure(WVCheckpointStatusCode::schemaMismatch,
                         "The output file omits required dynamic observer "
                         "state from its restart graph.",
                         file.destination);
      }
    }
    return validateObservationGraph();
  }

  WVCheckpointStatus validateDestinations(bool requireExisting) const {
    std::set<std::filesystem::path> observed;
    for (const auto &record : descriptor.record().outputFiles) {
      const auto destination =
          std::filesystem::absolute(record.destination).lexically_normal();
      if (!observed.insert(destination).second)
        return failure(WVCheckpointStatusCode::appendConflict,
                       "Output destinations are not unique.",
                       destination.string());
      std::error_code error;
      const bool exists = std::filesystem::exists(destination, error);
      if (error)
        return failure(WVCheckpointStatusCode::openFailure,
                       "Unable to inspect output destination: " +
                           error.message(),
                       destination.string());
      if (requireExisting != exists)
        return failure(requireExisting ? WVCheckpointStatusCode::openFailure
                                       : WVCheckpointStatusCode::commitFailure,
                       requireExisting
                           ? "Append destination does not exist."
                           : "New output destination already exists.",
                       destination.string());
      if (!std::filesystem::exists(destination.parent_path()))
        return failure(WVCheckpointStatusCode::openFailure,
                       "Output destination parent does not exist.",
                       destination.parent_path().string());
    }
    return WVCheckpointStatus::ok();
  }

  WVCheckpointStatus defineObserverMetadata(int parent,
                                            const WVObserverRecord &record,
                                            std::size_t ordinal,
                                            const std::string &path,
                                            const WVObservationSchema *schema,
                                            int &group) {
    const std::string groupName =
        "observingSystems-" + std::to_string(ordinal + 1);
    auto result = detail::checkedNetCDF(
        nc_def_grp(parent, groupName.c_str(), &group),
        "Observer metadata-group definition", path + "/" + groupName);
    if (!result)
      return result;
    const std::string observerPath = path + "/" + groupName;
    if (schema == nullptr) {
      result = detail::putTextAttribute(group, NC_GLOBAL, "AnnotatedClass",
                                        record.typeIdentifier, observerPath);
      if (result)
        result = detail::putTextAttribute(group, NC_GLOBAL,
                                          "portableIdentifier",
                                          record.identifier, observerPath);
      if (result)
        result = detail::putTextAttribute(group, NC_GLOBAL, "name",
                                          record.name, observerPath);
      if (!result)
        return result;
      int variable = -1;
      return defineScalar(group, "absTolerance", variable, observerPath);
    }
    const bool hadDeclaredAttributes = !schema->metadata.attributes.empty();
    const auto hasAttribute = [&](const std::string &name) {
      return std::any_of(schema->metadata.attributes.begin(),
                         schema->metadata.attributes.end(),
                         [&](const auto &attribute) {
                           return attribute.name == name;
                         });
    };
    for (const auto &[name, value] :
         std::array<std::pair<const char *, const std::string *>, 2>{
             {{"AnnotatedClass", &record.typeIdentifier},
              {"portableIdentifier", &record.identifier}}}) {
      if (hasAttribute(name))
        continue;
      result = detail::putTextAttribute(group, NC_GLOBAL, name, *value,
                                        observerPath);
      if (!result)
        return result;
    }
    if (!hasAttribute("name") &&
        (!hadDeclaredAttributes || !schema->preservesLegacyEncoding)) {
      result = detail::putTextAttribute(group, NC_GLOBAL, "name", record.name,
                                        observerPath);
      if (!result)
        return result;
    }
    for (const auto &attribute : schema->metadata.attributes) {
      result = detail::putTextAttribute(group, NC_GLOBAL,
                                        attribute.name.c_str(),
                                        attribute.value, observerPath);
      if (!result)
        return result;
    }
    for (const auto &attribute : schema->metadata.stringListAttributes) {
      result = putStringListAttribute(group, attribute.name.c_str(),
                                      attribute.values, observerPath);
      if (!result)
        return result;
    }
    if (!schema->preservesLegacyEncoding) {
      result = detail::putTextAttribute(
          group, NC_GLOBAL, "portableObservationSchemaIdentifier",
          schema->identifier, observerPath);
      if (!result)
        return result;
      result = detail::checkedNetCDF(
          nc_put_att_uint(group, NC_GLOBAL,
                          "portableObservationSchemaVersion", NC_UINT, 1,
                          &schema->version),
          "Observation-schema version definition", observerPath);
      if (!result)
        return result;
      result = detail::checkedNetCDF(
          nc_put_att_uint(group, NC_GLOBAL, "portableObserverContractVersion",
                          NC_UINT, 1, &record.contractVersion),
          "Observer contract-version definition", observerPath);
      if (!result)
        return result;
      std::vector<std::uint8_t> encodedConfiguration;
      const auto encoded = encodePortableTypedRecord(
          record.configuration, encodedConfiguration);
      if (!encoded)
        return failure(WVCheckpointStatusCode::schemaMismatch,
                       encoded.message, observerPath);
      result = detail::putTextAttribute(
          group, NC_GLOBAL, "portableObserverConfiguration",
          hexEncode(encodedConfiguration), observerPath);
      if (!result)
        return result;
      std::vector<std::uint8_t> schemaManifest;
      const auto manifestStatus =
          encodeObservationSchemaManifest(*schema, schemaManifest);
      if (!manifestStatus)
        return failure(WVCheckpointStatusCode::schemaMismatch,
                       manifestStatus.message, observerPath);
      result = detail::putTextAttribute(
          group, NC_GLOBAL, "portableObservationSchemaManifest",
          hexEncode(schemaManifest), observerPath);
      if (!result)
        return result;
    }
    for (const auto &metadataVariable : schema->metadata.variables) {
      int variable = -1;
      switch (metadataVariable.value.scalarType) {
      case WVObservationScalarType::real64:
        result = defineScalar(group, metadataVariable.name, variable,
                              observerPath);
        break;
      case WVObservationScalarType::complex64:
        result = detail::defineComplexVariable(
            group, metadataVariable.name, {}, observerPath);
        break;
      case WVObservationScalarType::integer64:
        result = detail::checkedNetCDF(
            nc_def_var(group, metadataVariable.name.c_str(), NC_INT64, 0,
                       nullptr, &variable),
            "Observer metadata-variable definition",
            observerPath + "/" + metadataVariable.name);
        break;
      case WVObservationScalarType::boolean8:
        result = detail::checkedNetCDF(
            nc_def_var(group, metadataVariable.name.c_str(), NC_UBYTE, 0,
                       nullptr, &variable),
            "Observer metadata-variable definition",
            observerPath + "/" + metadataVariable.name);
        if (result && metadataVariable.isLogicalType)
          result = detail::putByteAttribute(
              group, variable, "isLogicalType", 1,
              observerPath + "/" + metadataVariable.name);
        break;
      case WVObservationScalarType::text:
        result = detail::checkedNetCDF(
            nc_def_var(group, metadataVariable.name.c_str(), NC_STRING, 0,
                       nullptr, &variable),
            "Observer metadata-variable definition",
            observerPath + "/" + metadataVariable.name);
        break;
      }
      if (!result)
        return result;
    }
    return WVCheckpointStatus::ok();
  }

  WVCheckpointStatus defineGroup(File &file, Group &group,
                                 const std::array<int, 5> &rootDimensions) {
    auto result = detail::checkedNetCDF(
        nc_def_grp(file.id, group.record.name.c_str(), &group.id),
        "Output-group definition", "/" + group.record.name);
    if (!result)
      return result;
    const std::string path = "/" + group.record.name;
    result = detail::putTextAttribute(group.id, NC_GLOBAL, "AnnotatedClass",
                                      "WVModelOutputGroupEvenlySpaced", path);
    if (!result)
      return result;
    result = detail::putTextAttribute(group.id, NC_GLOBAL, "name",
                                      group.record.name, path);
    if (!result)
      return result;
    result = detail::putTextAttribute(group.id, NC_GLOBAL, "portableIdentifier",
                                      group.record.identifier, path);
    if (!result)
      return result;
    for (const char *name : {"outputInterval", "initialTime", "finalTime"}) {
      int variable = -1;
      result = defineScalar(group.id, name, variable, path);
      if (!result)
        return result;
    }
    if (!group.record.schedule.typeIdentifier.empty()) {
      result = detail::putTextAttribute(
          group.id, NC_GLOBAL, "portableScheduleTypeIdentifier",
          group.record.schedule.typeIdentifier, path);
      if (!result)
        return result;
      result = detail::checkedNetCDF(
          nc_put_att_uint(group.id, NC_GLOBAL,
                          "portableScheduleContractVersion", NC_UINT, 1,
                          &group.record.schedule.contractVersion),
          "Schedule-version attribute definition", path);
      if (!result)
        return result;
      std::vector<std::uint8_t> configuration;
      const auto encoded = encodePortableTypedRecord(
          group.record.schedule.configuration, configuration);
      if (!encoded)
        return failure(WVCheckpointStatusCode::invalidValue, encoded.message,
                       path);
      result = detail::putTextAttribute(group.id, NC_GLOBAL,
                                        "portableScheduleConfiguration",
                                        hexEncode(configuration), path);
      if (!result)
        return result;
    }
    int timeDimension = -1;
    result = detail::checkedNetCDF(
        nc_def_dim(group.id, "t", NC_UNLIMITED, &timeDimension),
        "Time-dimension definition", path + "/t");
    if (!result)
      return result;
    result = detail::defineDoubleVariable(group.id, "t", {timeDimension},
                                          group.timeId, path);
    if (!result)
      return result;
    for (const auto &attribute :
         std::array<std::pair<const char *, const char *>, 5>{
             {{"axis", "T"},
              {"calendar", "standard"},
              {"standard_name", "time"},
              {"long_name", "time"},
              {"units", "seconds since 1970-01-01 00:00:00"}}}) {
      result = detail::putTextAttribute(group.id, group.timeId, attribute.first,
                                        attribute.second, path + "/t");
      if (!result)
        return result;
    }
    if (!group.record.schedule.typeIdentifier.empty()) {
      result = detail::checkedNetCDF(
          nc_def_var(group.id, "portableScheduleOrdinal", NC_INT64, 1,
                     &timeDimension, &group.scheduleOrdinalId),
          "Schedule-ordinal variable definition", path);
      if (!result)
        return result;
      result = detail::checkedNetCDF(
          nc_def_var(group.id, "portableScheduleCursor", NC_STRING, 1,
                     &timeDimension, &group.scheduleCursorId),
          "Schedule-cursor variable definition", path);
      if (!result)
        return result;
    }
    if (group.record.containsCompleteCoefficientRestart) {
      group.coefficientReal.assign(stateLayout.coefficientFamilyCount(), -1);
      group.coefficientImag.assign(stateLayout.coefficientFamilyCount(), -1);
      for (std::size_t family = 0;
           family < stateLayout.coefficientFamilyCount(); ++family) {
        const auto &familyLayout = stateLayout.coefficientFamilies()[family];
        const std::string &name = familyLayout.identifier;
        std::vector<int> coefficientDimensions;
        if (!isDynamicsLinear)
          coefficientDimensions.push_back(timeDimension);
        if (familyLayout.spectralDimensions.size() == 2) {
          coefficientDimensions.push_back(rootDimensions[1]);
          coefficientDimensions.push_back(rootDimensions[0]);
        } else if (familyLayout.spectralDimensions.size() == 1) {
          coefficientDimensions.push_back(rootDimensions[1]);
        } else {
          return failure(WVCheckpointStatusCode::schemaMismatch,
                         "The transform coefficient family rank is not "
                         "supported by MATLAB NetCDF output.",
                         path + "/" + name);
        }
        result = detail::defineComplexVariable(group.id, name,
                                               coefficientDimensions, path);
        if (!result)
          return result;
        int real = -1;
        int imag = -1;
        result = detail::checkedNetCDF(
            nc_inq_varid(group.id, (name + "_real").c_str(), &real),
            "Coefficient-variable lookup", path + "/" + name + "_real");
        if (result)
          result = detail::checkedNetCDF(
              nc_inq_varid(group.id, (name + "_imag").c_str(), &imag),
              "Coefficient-variable lookup", path + "/" + name + "_imag");
        const auto *metadata = coefficientMetadata(name);
        if (metadata == nullptr)
          return failure(WVCheckpointStatusCode::schemaMismatch,
                         "Canonical coefficient metadata is unavailable.",
                         path + "/" + name);
        for (const int variable : {real, imag}) {
          if (result)
            result =
                detail::putTextAttribute(group.id, variable, "units",
                                         metadata->units, path + "/" + name);
          if (result)
            result = detail::putTextAttribute(group.id, variable, "long_name",
                                              metadata->description,
                                              path + "/" + name);
        }
        if (!result)
          return result;
        group.coefficientReal[family] = real;
        group.coefficientImag[family] = imag;
      }
    }

    int metadataRoot = -1;
    result = detail::checkedNetCDF(
        nc_def_grp(group.id, "observingSystems", &metadataRoot),
        "Observer metadata-root definition", path + "/observingSystems");
    if (!result)
      return result;

    std::map<std::string, WVObservationVariable> definedDerivedVariables;
    for (std::size_t ordinal = 0;
         ordinal < group.record.observerIdentifiers.size(); ++ordinal) {
      const auto *record = observer(group.record.observerIdentifiers[ordinal]);
      if (record == nullptr)
        return failure(WVCheckpointStatusCode::schemaMismatch,
                       "Output group references an unknown observer.", path);
      WVObservationSchema observationSchemaValue;
      const WVObservationSchema *observationSchema = nullptr;
      if (sampleSource != nullptr) {
        const auto sourceStatus = sampleSource->observationSchema(
            *record, observationSchemaValue);
        if (!sourceStatus)
          return failure(WVCheckpointStatusCode::unsupportedObserver,
                         sourceStatus.message, path);
        const auto schemaStatus =
            validateObservationSchema(observationSchemaValue);
        if (!schemaStatus)
          return failure(WVCheckpointStatusCode::schemaMismatch,
                         schemaStatus.message, path);
        observationSchema = &observationSchemaValue;
      }
      int observerMetadataGroup = -1;
      result = defineObserverMetadata(metadataRoot, *record, ordinal,
                                      path + "/observingSystems",
                                      observationSchema,
                                      observerMetadataGroup);
      if (!result)
        return result;

      if (observationSchema == nullptr)
        continue;
      for (const auto &specification : observationSchema->variables) {
        const bool canonicalCoefficient =
            group.record.containsCompleteCoefficientRestart &&
            coefficientMetadata(specification.name) != nullptr;
        if (canonicalCoefficient)
          continue;
        const auto existing = definedDerivedVariables.find(specification.name);
        if (existing != definedDerivedVariables.end()) {
          return failure(
              WVCheckpointStatusCode::schemaMismatch,
              "Observers declare duplicate persisted variable names without "
              "an alias contract.",
              path + "/" + specification.name);
        }
        definedDerivedVariables.emplace(specification.name, specification);
        std::vector<int> dimensions;
        if (specification.layout == WVObservationValueLayout::record)
          dimensions.push_back(timeDimension);
        std::vector<WVObservationAxis> variableAxes;
        variableAxes.reserve(specification.dimensionIdentifiers.size());
        for (const auto &identifier : specification.dimensionIdentifiers) {
          const auto *axisDefinition = axis(*observationSchema, identifier);
          if (axisDefinition == nullptr)
            return failure(WVCheckpointStatusCode::schemaMismatch,
                           "Observation variable references an unknown axis.",
                           path + "/" + specification.name);
          variableAxes.push_back(*axisDefinition);
        }
        for (std::size_t reverse = variableAxes.size(); reverse > 0; --reverse) {
          const auto logicalIndex = reverse - 1;
          int dimension = -1;
          const auto &axisDefinition = variableAxes[logicalIndex];
          const auto &name = axisDefinition.name;
          if (nc_inq_dimid(group.id, name.c_str(), &dimension) != NC_NOERR) {
            result = detail::checkedNetCDF(
                nc_def_dim(group.id, name.c_str(),
                           axisDefinition.kind == WVObservationAxisKind::unlimited
                               ? NC_UNLIMITED
                               : axisDefinition.extent,
                           &dimension),
                "Observer dimension definition", path + "/" + name);
            if (!result)
              return result;
          }
          if (axisDefinition.kind == WVObservationAxisKind::unlimited) {
            const auto foundAxis = std::find_if(
                group.unlimitedAxes.begin(), group.unlimitedAxes.end(),
                [&](const auto &candidate) { return candidate.name == name; });
            if (foundAxis == group.unlimitedAxes.end()) {
              UnlimitedAxis unlimited;
              unlimited.name = name;
              unlimited.dimensionId = dimension;
              const std::string progressName = "portableCommitted_" + name;
              result = detail::checkedNetCDF(
                  nc_def_var(group.id, progressName.c_str(), NC_INT64, 1,
                             &timeDimension, &unlimited.progressVariableId),
                  "Observation-axis progress definition",
                  path + "/" + progressName);
              if (!result)
                return result;
              group.unlimitedAxes.push_back(std::move(unlimited));
            } else if (foundAxis->dimensionId != dimension) {
              return failure(WVCheckpointStatusCode::schemaMismatch,
                             "Observation schemas disagree on an unlimited axis.",
                             path + "/" + name);
            }
          }
          dimensions.push_back(dimension);
        }
        Variable variable;
        variable.observerIdentifier = record->identifier;
        variable.schemaIdentifier = observationSchema->identifier;
        variable.schemaVersion = observationSchema->version;
        variable.specification = specification;
        variable.axes = std::move(variableAxes);
        if (specification.scalarType == WVObservationScalarType::complex64) {
          result = detail::defineComplexVariable(group.id, specification.name,
                                                 dimensions, path);
          if (!result)
            return result;
        } else if (specification.scalarType ==
                   WVObservationScalarType::real64) {
          result = detail::defineDoubleVariable(
              group.id, specification.name, dimensions, variable.realId, path);
          if (!result)
            return result;
        } else {
          const nc_type type =
              specification.scalarType == WVObservationScalarType::integer64
                  ? NC_INT64
                  : specification.scalarType ==
                            WVObservationScalarType::boolean8
                        ? NC_UBYTE
                        : NC_STRING;
          result = detail::checkedNetCDF(
              nc_def_var(group.id, specification.name.c_str(), type,
                         static_cast<int>(dimensions.size()),
                         dimensions.empty() ? nullptr : dimensions.data(),
                         &variable.realId),
              "Observer-variable definition",
              path + "/" + specification.name);
          if (!result)
            return result;
          if (specification.scalarType == WVObservationScalarType::boolean8) {
            result = detail::putByteAttribute(
                group.id, variable.realId, "isLogicalType", 1,
                path + "/" + specification.name);
            if (!result)
              return result;
          }
        }
        const auto applyAttributes = [&](int variableId,
                                         const std::string &variablePath) {
          auto attributeResult = WVCheckpointStatus::ok();
          if (!specification.units.empty())
            attributeResult =
                detail::putTextAttribute(group.id, variableId, "units",
                                         specification.units, variablePath);
          if (attributeResult && !specification.description.empty())
            attributeResult =
                detail::putTextAttribute(group.id, variableId, "long_name",
                                         specification.description, variablePath);
          for (const auto &attribute : specification.attributes) {
            if (!attributeResult)
              break;
            attributeResult = detail::putTextAttribute(
                group.id, variableId, attribute.name.c_str(), attribute.value,
                variablePath);
          }
          if (attributeResult && !observationSchema->preservesLegacyEncoding) {
            attributeResult = detail::putTextAttribute(
                group.id, variableId, "portableObservationVariableIdentifier",
                specification.identifier, variablePath);
          }
          if (attributeResult && !observationSchema->preservesLegacyEncoding)
            attributeResult = detail::putTextAttribute(
                group.id, variableId,
                "portableObservationObserverIdentifier", record->identifier,
                variablePath);
          if (attributeResult && !observationSchema->preservesLegacyEncoding)
            attributeResult = detail::putTextAttribute(
                group.id, variableId, "portableObservationSchemaIdentifier",
                observationSchema->identifier, variablePath);
          if (attributeResult && !observationSchema->preservesLegacyEncoding)
            attributeResult = detail::checkedNetCDF(
                nc_put_att_uint(group.id, variableId,
                                "portableObservationSchemaVersion", NC_UINT,
                                1, &observationSchema->version),
                "Observation-schema version definition", variablePath);
          if (attributeResult && !observationSchema->preservesLegacyEncoding)
            attributeResult = detail::putTextAttribute(
                group.id, variableId, "portableObservationValueLayout",
                valueLayoutName(specification.layout), variablePath);
          if (attributeResult && !observationSchema->preservesLegacyEncoding &&
              !specification.dimensionIdentifiers.empty()) {
            std::vector<const char *> rawDimensions;
            rawDimensions.reserve(specification.dimensionIdentifiers.size());
            for (const auto &identifier :
                 specification.dimensionIdentifiers)
              rawDimensions.push_back(identifier.c_str());
            attributeResult = detail::checkedNetCDF(
                nc_put_att_string(
                    group.id, variableId,
                    "portableObservationDimensionIdentifiers",
                    rawDimensions.size(), rawDimensions.data()),
                "Observation dimension-identifier definition", variablePath);
            if (attributeResult) {
              std::vector<const char *> rawRoles;
              std::vector<std::string> roles;
              roles.reserve(variable.axes.size());
              rawRoles.reserve(variable.axes.size());
              for (const auto &axisDefinition : variable.axes)
                roles.emplace_back(
                    coordinateRoleName(axisDefinition.coordinateRole));
              for (const auto &role : roles)
                rawRoles.push_back(role.c_str());
              attributeResult = detail::checkedNetCDF(
                  nc_put_att_string(
                      group.id, variableId,
                      "portableObservationAxisCoordinateRoles",
                      rawRoles.size(), rawRoles.data()),
                  "Observation axis-role definition", variablePath);
            }
          }
          if (attributeResult && !observationSchema->preservesLegacyEncoding &&
              specification.coordinateRole !=
                  WVObservationCoordinateRole::none)
            attributeResult = detail::putTextAttribute(
                group.id, variableId, "portableCoordinateRole",
                coordinateRoleName(specification.coordinateRole), variablePath);
          if (attributeResult && !observationSchema->preservesLegacyEncoding &&
              specification.raggedRole != WVObservationRaggedRole::none) {
            attributeResult = detail::putTextAttribute(
                group.id, variableId, "portableRaggedRole",
                raggedRoleName(specification.raggedRole), variablePath);
            if (attributeResult)
              attributeResult = detail::putTextAttribute(
                  group.id, variableId, "portableRaggedChildAxis",
                  specification.raggedChildAxisIdentifier, variablePath);
          }
          return attributeResult;
        };
        if (specification.scalarType == WVObservationScalarType::complex64) {
          int real = -1;
          int imag = -1;
          result = detail::checkedNetCDF(
              nc_inq_varid(group.id, (specification.name + "_real").c_str(),
                           &real),
              "Observer-variable lookup", path);
          if (result)
            result = detail::checkedNetCDF(
                nc_inq_varid(group.id, (specification.name + "_imag").c_str(),
                             &imag),
                "Observer-variable lookup", path);
          variable.realId = real;
          variable.imagId = imag;
          if (result)
            result = applyAttributes(real,
                                     path + "/" + specification.name + "_real");
          if (result)
            result = applyAttributes(imag,
                                     path + "/" + specification.name + "_imag");
        } else {
          result =
              applyAttributes(variable.realId, path + "/" + specification.name);
        }
        if (!result)
          return result;
        group.derivedVariables.push_back(std::move(variable));
      }
      group.observerSchemas.push_back(
          {record->identifier, std::move(observationSchemaValue)});
    }
    const auto bindingStatus = compileVariableWriteBindings(group);
    if (!bindingStatus)
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     bindingStatus.message, path);
    return WVCheckpointStatus::ok();
  }

  WVCheckpointStatus writeObserverMetadata(int fileId, const File &file) {
    for (const auto &group : file.groups) {
      const std::string path = "/" + group.record.name;
      auto result = writeScalar(group.id, "outputInterval",
                                group.record.schedule.outputInterval, path);
      if (!result)
        return result;
      result = writeScalar(group.id, "initialTime",
                           group.record.schedule.initialTime, path);
      if (!result)
        return result;
      result = writeScalar(group.id, "finalTime",
                           group.record.schedule.finalTime, path);
      if (!result)
        return result;
      int metadataRoot = -1;
      result = detail::checkedNetCDF(
          nc_inq_ncid(group.id, "observingSystems", &metadataRoot),
          "Observer metadata-root lookup", path + "/observingSystems");
      if (!result)
        return result;
      for (std::size_t ordinal = 0;
           ordinal < group.record.observerIdentifiers.size(); ++ordinal) {
        const auto *record =
            observer(group.record.observerIdentifiers[ordinal]);
        const std::string groupName =
            "observingSystems-" + std::to_string(ordinal + 1);
        int metadata = -1;
        result = detail::checkedNetCDF(
            nc_inq_ncid(metadataRoot, groupName.c_str(), &metadata),
            "Observer metadata-group lookup",
            path + "/observingSystems/" + groupName);
        if (!result)
          return result;
        const std::string metadataPath =
            path + "/observingSystems/" + groupName;
        const auto observerSchema = std::find_if(
            group.observerSchemas.begin(), group.observerSchemas.end(),
            [&](const auto &candidate) {
              return candidate.observerIdentifier == record->identifier;
            });
        if (observerSchema == group.observerSchemas.end()) {
          const auto *block = stateLayout.coefficientFamilies().empty()
                                  ? nullptr
                                  : stateBlockRecord(
                                        stateLayout.coefficientFamilies()
                                            .front()
                                            .identifier);
          result = writeScalar(
              metadata, "absTolerance",
              block == nullptr ? 1e-6 : block->absoluteTolerance, metadataPath);
          if (!result)
            return result;
          continue;
        }
        for (const auto &metadataVariable :
             observerSchema->schema.metadata.variables) {
          const auto &value = metadataVariable.value;
          if (value.elementCount() != 1)
            return failure(WVCheckpointStatusCode::shapeMismatch,
                           "Observer metadata variables must be scalar.",
                           metadataPath + "/" + metadataVariable.name);
          int variable = -1;
          if (value.scalarType == WVObservationScalarType::complex64) {
            int real = -1;
            int imaginary = -1;
            result = detail::checkedNetCDF(
                nc_inq_varid(metadata,
                             (metadataVariable.name + "_real").c_str(), &real),
                "Observer metadata-variable lookup", metadataPath);
            if (result)
              result = detail::checkedNetCDF(
                  nc_inq_varid(metadata,
                               (metadataVariable.name + "_imag").c_str(),
                               &imaginary),
                  "Observer metadata-variable lookup", metadataPath);
            if (result) {
              const auto entry = value.complex64Data()[0];
              result = detail::checkedNetCDF(
                  nc_put_var_double(metadata, real, &entry.real),
                  "Observer metadata-variable write", metadataPath);
              if (result)
                result = detail::checkedNetCDF(
                    nc_put_var_double(metadata, imaginary, &entry.imag),
                    "Observer metadata-variable write", metadataPath);
            }
          } else {
            result = detail::checkedNetCDF(
                nc_inq_varid(metadata, metadataVariable.name.c_str(),
                             &variable),
                "Observer metadata-variable lookup",
                metadataPath + "/" + metadataVariable.name);
            if (result &&
                value.scalarType == WVObservationScalarType::real64)
              result = detail::checkedNetCDF(
                  nc_put_var_double(metadata, variable, value.real64Data()),
                  "Observer metadata-variable write", metadataPath);
            else if (result &&
                     value.scalarType == WVObservationScalarType::integer64) {
              const long long entry =
                  static_cast<long long>(value.integer64Data()[0]);
              result = detail::checkedNetCDF(
                  nc_put_var_longlong(metadata, variable, &entry),
                  "Observer metadata-variable write", metadataPath);
            } else if (result &&
                       value.scalarType == WVObservationScalarType::boolean8)
              result = detail::checkedNetCDF(
                  nc_put_var_uchar(metadata, variable, value.boolean8Data()),
                  "Observer metadata-variable write", metadataPath);
            else if (result &&
                     value.scalarType == WVObservationScalarType::text) {
              const char *entry = value.textData()[0].c_str();
              result = detail::checkedNetCDF(
                  nc_put_var_string(metadata, variable, &entry),
                  "Observer metadata-variable write", metadataPath);
            }
          }
          if (!result)
            return result;
        }
      }
    }
    (void)fileId;
    return WVCheckpointStatus::ok();
  }

  WVCheckpointStatus writeObservationValue(
      Group &group, Variable &variable,
      const WVObservationValue &value, std::size_t recordIndex,
      const std::string &path) {
    auto &start = variable.writeStart;
    auto &count = variable.writeCount;
    const auto expectedRank =
        variable.axes.size() +
        (variable.specification.layout == WVObservationValueLayout::record
             ? 1
             : 0);
    if (start.size() != expectedRank || count.size() != expectedRank ||
        variable.unlimitedAxisIndices.size() != variable.axes.size())
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     "Observation write bindings were not compiled.", path);
    std::fill(start.begin(), start.end(), 0);
    std::fill(count.begin(), count.end(), 1);
    std::size_t persistedDimension = 0;
    if (variable.specification.layout == WVObservationValueLayout::record) {
      start[persistedDimension] = recordIndex;
      ++persistedDimension;
    }
    for (std::size_t reverse = variable.axes.size(); reverse > 0; --reverse) {
      const auto logicalIndex = reverse - 1;
      const auto &axisDefinition = variable.axes[logicalIndex];
      std::size_t axisStart = 0;
      if (axisDefinition.kind == WVObservationAxisKind::unlimited) {
        const auto axisIndex = variable.unlimitedAxisIndices[logicalIndex];
        if (axisIndex >= group.unlimitedAxes.size())
          return failure(WVCheckpointStatusCode::schemaMismatch,
                         "Observation batch references an unresolved unlimited axis.",
                         path);
        axisStart = group.unlimitedAxes[axisIndex].committedCount;
      }
      start[persistedDimension] = axisStart;
      count[persistedDimension] = value.extents[logicalIndex];
      ++persistedDimension;
    }
    const auto elementCount = value.elementCount();
    if (elementCount == 0)
      return WVCheckpointStatus::ok();
    const auto putReal = [&](int variableId, const double *values,
                             const std::string &variablePath) {
      return detail::checkedNetCDF(
          start.empty()
              ? nc_put_var_double(group.id, variableId, values)
              : nc_put_vara_double(group.id, variableId, start.data(),
                                   count.data(), values),
          "Observation value write", variablePath);
    };
    if (value.scalarType == WVObservationScalarType::real64)
      return putReal(variable.realId, value.real64Data(), path);
    if (value.scalarType == WVObservationScalarType::complex64) {
      group.realWriteScratch.resize(elementCount);
      group.imaginaryWriteScratch.resize(elementCount);
      for (std::size_t index = 0; index < elementCount; ++index) {
        group.realWriteScratch[index] = value.complex64Data()[index].real;
        group.imaginaryWriteScratch[index] =
            value.complex64Data()[index].imag;
      }
      auto result = putReal(variable.realId, group.realWriteScratch.data(),
                            path + "_real");
      if (!result)
        return result;
      return putReal(variable.imagId, group.imaginaryWriteScratch.data(),
                     path + "_imag");
    }
    if (value.scalarType == WVObservationScalarType::integer64) {
      group.integerWriteScratch.resize(elementCount);
      std::transform(value.integer64Data(),
                     value.integer64Data() + elementCount,
                     group.integerWriteScratch.begin(),
                     [](std::int64_t entry) {
                       return static_cast<long long>(entry);
                     });
      return detail::checkedNetCDF(
          start.empty()
              ? nc_put_var_longlong(group.id, variable.realId,
                                    group.integerWriteScratch.data())
              : nc_put_vara_longlong(group.id, variable.realId, start.data(),
                                     count.data(),
                                     group.integerWriteScratch.data()),
          "Integer observation value write", path);
    }
    if (value.scalarType == WVObservationScalarType::boolean8)
      return detail::checkedNetCDF(
          start.empty()
              ? nc_put_var_uchar(group.id, variable.realId,
                                 value.boolean8Data())
              : nc_put_vara_uchar(group.id, variable.realId, start.data(),
                                  count.data(), value.boolean8Data()),
          "Boolean observation value write", path);
    group.textWriteScratch.resize(elementCount);
    for (std::size_t index = 0; index < elementCount; ++index)
      group.textWriteScratch[index] = value.textData()[index].c_str();
    return detail::checkedNetCDF(
        start.empty()
            ? nc_put_var_string(group.id, variable.realId,
                                group.textWriteScratch.data())
            : nc_put_vara_string(group.id, variable.realId, start.data(),
                                 count.data(), group.textWriteScratch.data()),
        "Text observation value write", path);
  }

  WVCheckpointStatus writeInitialAndStaticValues(File &file) {
    if (constructionCheckpoint == nullptr)
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     "Initial output has no checkpoint template.", "/");
    for (auto &group : file.groups) {
      const std::string path = "/" + group.record.name;
      if (group.record.containsCompleteCoefficientRestart &&
          isDynamicsLinear) {
        for (std::size_t family = 0;
             family < stateLayout.coefficientFamilyCount(); ++family) {
          const auto &familyLayout = stateLayout.coefficientFamilies()[family];
          const std::string &name = familyLayout.identifier;
          const auto coefficientCount = familyLayout.elementCount;
          const auto *values = constructionCoefficientData(family);
          if (values == nullptr)
            return failure(WVCheckpointStatusCode::shapeMismatch,
                           "Initial coefficient storage is incomplete.",
                           path + "/" + name);
          int real = -1;
          int imag = -1;
          auto result = detail::checkedNetCDF(
              nc_inq_varid(group.id, (name + "_real").c_str(), &real),
              "Initial coefficient-variable lookup", path + "/" + name);
          if (result)
            result = detail::checkedNetCDF(
                nc_inq_varid(group.id, (name + "_imag").c_str(), &imag),
                "Initial coefficient-variable lookup", path + "/" + name);
          if (!result)
            return result;
          std::vector<double> realValues(coefficientCount);
          std::vector<double> imaginaryValues(coefficientCount);
          for (std::size_t index = 0; index < coefficientCount; ++index) {
            realValues[index] = values[index].real;
            imaginaryValues[index] = values[index].imag;
          }
          result = detail::checkedNetCDF(
              nc_put_var_double(group.id, real, realValues.data()),
              "Initial coefficient write", path + "/" + name + "_real");
          if (result)
            result = detail::checkedNetCDF(
                nc_put_var_double(group.id, imag, imaginaryValues.data()),
                "Initial coefficient write", path + "/" + name + "_imag");
          if (!result)
            return result;
        }
      }
      for (const auto &observerSchema : group.observerSchemas) {
        const auto *record = observer(observerSchema.observerIdentifier);
        WVObservationBatch batch;
        const auto sourceStatus =
            sampleSource->initialObservationBatch(*record, batch);
        if (!sourceStatus)
          return failure(WVCheckpointStatusCode::unsupportedObserver,
                         sourceStatus.message, path);
        const auto batchStatus =
            validateObservationBatch(observerSchema.schema, batch);
        if (!batchStatus)
          return failure(WVCheckpointStatusCode::shapeMismatch,
                         batchStatus.message, path);
        const auto batchMetrics = batch.metrics();
        metrics.batchRetainedStorageBytes = std::max(
            metrics.batchRetainedStorageBytes,
            batchMetrics.retainedStorageBytes);
        metrics.batchMaximumLiveBytes =
            std::max(metrics.batchMaximumLiveBytes, batchMetrics.liveBytes);
        for (const auto &value : batch.values) {
          const auto variable = std::find_if(
              group.derivedVariables.begin(), group.derivedVariables.end(),
              [&](const auto &candidate) {
                return candidate.observerIdentifier ==
                           observerSchema.observerIdentifier &&
                       candidate.specification.identifier ==
                           value.variableIdentifier;
              });
          if (variable == group.derivedVariables.end()) {
            const auto declared = std::find_if(
                observerSchema.schema.variables.begin(),
                observerSchema.schema.variables.end(), [&](const auto &item) {
                  return item.identifier == value.variableIdentifier;
                });
            if (declared != observerSchema.schema.variables.end() &&
                group.record.containsCompleteCoefficientRestart &&
                coefficientMetadata(declared->name) != nullptr)
              continue;
            return failure(WVCheckpointStatusCode::schemaMismatch,
                           "Initial observation value has no persistence variable.",
                           path);
          }
          const auto result = writeObservationValue(
              group, *variable, value, 0,
              path + "/" + variable->specification.name);
          if (!result)
            return result;
        }
      }
    }
    return WVCheckpointStatus::ok();
  }

  WVCheckpointStatus stageFiles() {
    if (sampleSource != nullptr) {
      std::vector<WVCoefficientFamilyConstView> coefficientViews;
      coefficientViews.reserve(stateLayout.coefficientFamilyCount());
      for (std::size_t family = 0;
           family < stateLayout.coefficientFamilyCount(); ++family)
        coefficientViews.push_back(
            {&stateLayout.coefficientFamilies()[family],
             constructionCoefficientData(family)});
      WVIntegrationState initial;
      transformAdapter->bindConstructionState(*constructionCheckpoint,
                                              initial);
      initial.coefficientFamilies = coefficientViews.data();
      initial.coefficientFamilyCount = coefficientViews.size();
      const auto status = sampleSource->prepareInitial(initial);
      if (!status)
        return failure(WVCheckpointStatusCode::unsupportedObserver,
                       status.message, "/observingSystems");
    }
    files.clear();
    files.reserve(descriptor.record().outputFiles.size());
    for (const auto &record : descriptor.record().outputFiles) {
      File file(record);
      file.destination =
          std::filesystem::absolute(record.destination).lexically_normal();
      file.staging = temporaryPath(file.destination);
      int id = -1;
      auto result = detail::checkedNetCDF(
          nc_create(file.staging.c_str(), NC_NETCDF4 | NC_NOCLOBBER, &id),
          "Output-file creation", file.staging.string());
      if (!result)
        return result;
      file.id = id;
      files.push_back(std::move(file));
      auto &staged = files.back();
      std::array<int, 5> rootDimensions{};
      std::vector<int> forcingGroups;
      std::vector<const WVFrozenForcingEntry *> forcingEntries;
      result = transformAdapter->defineRoot(
          staged.id, *constructionCheckpoint, catalog->forcings(),
          isDynamicsLinear, rootDimensions, forcingGroups,
          forcingEntries);
      if (!result)
        return result;
      result = detail::putTextAttribute(staged.id, NC_GLOBAL,
                                        "portableFileIdentifier",
                                        record.identifier, "/");
      if (!result)
        return result;
      staged.groups.reserve(record.groups.size());
      for (const auto &groupRecord : record.groups) {
        Group group(groupRecord);
        result = defineGroup(staged, group, rootDimensions);
        if (!result)
          return result;
        staged.groups.push_back(std::move(group));
      }
      result = detail::checkedNetCDF(nc_enddef(staged.id),
                                     "End output-file definition",
                                     staged.staging.string());
      if (!result)
        return result;
      result = transformAdapter->writeRoot(
          staged.id, *constructionCheckpoint, catalog->forcings(),
          forcingGroups, forcingEntries);
      if (!result)
        return result;
      result = writeObserverMetadata(staged.id, staged);
      if (!result)
        return result;
      result = writeInitialAndStaticValues(staged);
      if (!result)
        return result;
      result = detail::checkedNetCDF(nc_sync(staged.id),
                                     "Output-file initialization sync",
                                     staged.staging.string());
      if (!result)
        return result;
    }
    return WVCheckpointStatus::ok();
  }

  WVCheckpointStatus commitStagedFiles() {
    for (auto &file : files) {
      if (file.id >= 0) {
        const int id = std::exchange(file.id, -1);
        const auto result = detail::checkedNetCDF(
            nc_close(id), "Output-file initialization close",
            file.staging.string());
        if (!result)
          return result;
      }
    }
    std::vector<std::filesystem::path> committed;
    for (auto &file : files) {
      std::error_code error;
      std::filesystem::create_hard_link(file.staging, file.destination, error);
      if (error) {
        for (const auto &path : committed) {
          std::error_code ignored;
          std::filesystem::remove(path, ignored);
        }
        return failure(WVCheckpointStatusCode::commitFailure,
                       "Unable to commit output file set: " + error.message(),
                       file.destination.string());
      }
      committed.push_back(file.destination);
    }
    for (auto &file : files) {
      std::error_code ignored;
      std::filesystem::remove(file.staging, ignored);
      int id = -1;
      auto result = detail::checkedNetCDF(
          nc_open(file.destination.c_str(), NC_WRITE, &id),
          "Committed output-file open", file.destination.string());
      if (!result) {
        for (auto &candidate : files) {
          if (candidate.id >= 0) {
            nc_close(std::exchange(candidate.id, -1));
          }
        }
        for (const auto &path : committed) {
          std::error_code rollbackError;
          std::filesystem::remove(path, rollbackError);
        }
        return result;
      }
      file.id = id;
      result = resolveFile(file);
      if (!result) {
        for (auto &candidate : files) {
          if (candidate.id >= 0)
            nc_close(std::exchange(candidate.id, -1));
        }
        for (const auto &path : committed) {
          std::error_code rollbackError;
          std::filesystem::remove(path, rollbackError);
        }
        return result;
      }
    }
    return WVCheckpointStatus::ok();
  }

  WVCheckpointStatus replaceStagedFiles() {
    for (auto &file : files) {
      if (file.id >= 0) {
        const int id = std::exchange(file.id, -1);
        const auto result = detail::checkedNetCDF(
            nc_close(id), "Output-file initialization close",
            file.staging.string());
        if (!result)
          return result;
      }
    }

    struct Replacement {
      std::filesystem::path destination;
      std::filesystem::path backup;
      bool moved = false;
      bool installed = false;
    };
    std::vector<Replacement> replacements;
    replacements.reserve(files.size());
    for (const auto &file : files)
      replacements.push_back(
          {file.destination, temporaryPath(file.destination), false, false});

    auto rollback = [&]() noexcept {
      for (auto &file : files) {
        if (file.id >= 0)
          nc_close(std::exchange(file.id, -1));
      }
      for (auto iterator = replacements.rbegin();
           iterator != replacements.rend(); ++iterator) {
        std::error_code ignored;
        if (iterator->installed)
          std::filesystem::remove(iterator->destination, ignored);
        if (iterator->moved)
          std::filesystem::rename(iterator->backup, iterator->destination,
                                  ignored);
      }
    };

    for (auto &replacement : replacements) {
      std::error_code error;
      std::filesystem::rename(replacement.destination, replacement.backup,
                              error);
      if (error) {
        rollback();
        return failure(WVCheckpointStatusCode::commitFailure,
                       "Unable to stage the existing output set for "
                       "replacement: " +
                           error.message(),
                       replacement.destination.string());
      }
      replacement.moved = true;
    }
    for (std::size_t index = 0; index < files.size(); ++index) {
      std::error_code error;
      std::filesystem::create_hard_link(files[index].staging,
                                        files[index].destination, error);
      if (error) {
        rollback();
        return failure(WVCheckpointStatusCode::commitFailure,
                       "Unable to install the replacement output set: " +
                           error.message(),
                       files[index].destination.string());
      }
      replacements[index].installed = true;
    }

    for (auto &file : files) {
      int id = -1;
      auto result = detail::checkedNetCDF(
          nc_open(file.destination.c_str(), NC_WRITE, &id),
          "Replacement output-file open", file.destination.string());
      if (!result) {
        rollback();
        return result;
      }
      file.id = id;
      result = resolveFile(file);
      if (!result) {
        rollback();
        return result;
      }
    }

    // Replacement is committed once every new file has opened and resolved.
    // Backup cleanup is best-effort: a cleanup failure must not report a
    // rollback-safe failure after an earlier backup has already been removed.
    for (const auto &replacement : replacements) {
      std::error_code ignored;
      std::filesystem::remove(replacement.backup, ignored);
    }
    for (auto &file : files) {
      std::error_code ignored;
      std::filesystem::remove(file.staging, ignored);
    }
    return WVCheckpointStatus::ok();
  }

  WVCheckpointStatus resolveFile(File &file) {
    for (auto &group : file.groups) {
      auto result = detail::checkedNetCDF(
          nc_inq_ncid(file.id, group.record.name.c_str(), &group.id),
          "Output-group lookup", "/" + group.record.name);
      if (!result)
        return result;
      for (const auto &expected :
           std::array<std::pair<const char *, double>, 3>{
               {{"outputInterval", group.record.schedule.outputInterval},
                {"initialTime", group.record.schedule.initialTime},
                {"finalTime", group.record.schedule.finalTime}}}) {
        double observed = 0.0;
        result = detail::readDoubleScalar(group.id, expected.first, observed,
                                          "/" + group.record.name);
        if (!result)
          return result;
        if (observed != expected.second)
          return failure(WVCheckpointStatusCode::appendConflict,
                         "Existing output schedule does not match append "
                         "configuration.",
                         "/" + group.record.name + "/" + expected.first);
      }
      if (!group.record.schedule.typeIdentifier.empty()) {
        std::string observedType;
        result = detail::readTextAttribute(
            group.id, "portableScheduleTypeIdentifier", observedType,
            "/" + group.record.name);
        if (!result || observedType != group.record.schedule.typeIdentifier)
          return failure(WVCheckpointStatusCode::appendConflict,
                         "Existing algorithmic schedule identity differs.",
                         "/" + group.record.name);
        unsigned int observedVersion = 0;
        result = detail::checkedNetCDF(
            nc_get_att_uint(group.id, NC_GLOBAL,
                            "portableScheduleContractVersion",
                            &observedVersion),
            "Schedule-version attribute read", "/" + group.record.name);
        if (!result || observedVersion != group.record.schedule.contractVersion)
          return failure(WVCheckpointStatusCode::appendConflict,
                         "Existing algorithmic schedule version differs.",
                         "/" + group.record.name);
        std::vector<std::uint8_t> configuration;
        const auto encoded = encodePortableTypedRecord(
            group.record.schedule.configuration, configuration);
        if (!encoded)
          return failure(WVCheckpointStatusCode::invalidValue, encoded.message,
                         "/" + group.record.name);
        std::string observedConfiguration;
        result = detail::readTextAttribute(
            group.id, "portableScheduleConfiguration", observedConfiguration,
            "/" + group.record.name);
        if (!result || observedConfiguration != hexEncode(configuration))
          return failure(WVCheckpointStatusCode::appendConflict,
                         "Existing algorithmic schedule configuration differs.",
                         "/" + group.record.name);
      }
      if (!group.observerSchemas.empty()) {
        int metadataRoot = -1;
        result = detail::checkedNetCDF(
            nc_inq_ncid(group.id, "observingSystems", &metadataRoot),
            "Observer metadata-root lookup",
            "/" + group.record.name + "/observingSystems");
        if (!result)
          return result;
        for (std::size_t ordinal = 0;
             ordinal < group.record.observerIdentifiers.size(); ++ordinal) {
          const auto schema = std::find_if(
              group.observerSchemas.begin(), group.observerSchemas.end(),
              [&](const auto &candidate) {
                return candidate.observerIdentifier ==
                       group.record.observerIdentifiers[ordinal];
              });
          if (schema == group.observerSchemas.end() ||
              schema->schema.preservesLegacyEncoding)
            continue;
          int metadata = -1;
          const std::string metadataName =
              "observingSystems-" + std::to_string(ordinal + 1);
          result = detail::checkedNetCDF(
              nc_inq_ncid(metadataRoot, metadataName.c_str(), &metadata),
              "Observer metadata-group lookup",
              "/" + group.record.name + "/observingSystems/" + metadataName);
          if (!result)
            return result;
          std::string observedIdentifier;
          result = detail::readTextAttribute(
              metadata, "portableObservationSchemaIdentifier",
              observedIdentifier,
              "/" + group.record.name + "/observingSystems/" + metadataName);
          if (!result || observedIdentifier != schema->schema.identifier)
            return failure(WVCheckpointStatusCode::appendConflict,
                           "Observation schema identity changed.",
                           "/" + group.record.name + "/observingSystems/" +
                               metadataName);
          unsigned int observedVersion = 0;
          result = detail::checkedNetCDF(
              nc_get_att_uint(metadata, NC_GLOBAL,
                              "portableObservationSchemaVersion",
                              &observedVersion),
              "Observation-schema version read",
              "/" + group.record.name + "/observingSystems/" + metadataName);
          if (!result || observedVersion != schema->schema.version)
            return failure(WVCheckpointStatusCode::appendConflict,
                           "Observation schema version changed.",
                           "/" + group.record.name + "/observingSystems/" +
                               metadataName);
          const auto *observerRecord = observer(schema->observerIdentifier);
          if (observerRecord == nullptr)
            return failure(WVCheckpointStatusCode::schemaMismatch,
                           "Observation schema has no resolved observer.",
                           "/" + group.record.name + "/observingSystems/" +
                               metadataName);
          unsigned int observedObserverVersion = 0;
          result = detail::checkedNetCDF(
              nc_get_att_uint(metadata, NC_GLOBAL,
                              "portableObserverContractVersion",
                              &observedObserverVersion),
              "Observer contract-version read",
              "/" + group.record.name + "/observingSystems/" + metadataName);
          if (!result ||
              observedObserverVersion != observerRecord->contractVersion)
            return failure(WVCheckpointStatusCode::appendConflict,
                           "Observer contract version changed.",
                           "/" + group.record.name + "/observingSystems/" +
                               metadataName);
          std::vector<std::uint8_t> configurationBytes;
          const auto encoded = encodePortableTypedRecord(
              observerRecord->configuration, configurationBytes);
          if (!encoded)
            return failure(WVCheckpointStatusCode::schemaMismatch,
                           encoded.message, "/" + group.record.name);
          std::string observedConfiguration;
          result = detail::readTextAttribute(
              metadata, "portableObserverConfiguration",
              observedConfiguration,
              "/" + group.record.name + "/observingSystems/" + metadataName);
          if (!result || observedConfiguration != hexEncode(configurationBytes))
            return failure(WVCheckpointStatusCode::appendConflict,
                           "Observer typed configuration changed.",
                           "/" + group.record.name + "/observingSystems/" +
                               metadataName);
          std::vector<std::uint8_t> manifestBytes;
          const auto manifestStatus = encodeObservationSchemaManifest(
              schema->schema, manifestBytes);
          if (!manifestStatus)
            return failure(WVCheckpointStatusCode::schemaMismatch,
                           manifestStatus.message,
                           "/" + group.record.name);
          std::string observedManifest;
          result = detail::readTextAttribute(
              metadata, "portableObservationSchemaManifest",
              observedManifest,
              "/" + group.record.name + "/observingSystems/" + metadataName);
          if (!result || observedManifest != hexEncode(manifestBytes))
            return failure(WVCheckpointStatusCode::appendConflict,
                           "Observation schema manifest changed.",
                           "/" + group.record.name + "/observingSystems/" +
                               metadataName);
        }
      }
      result = detail::checkedNetCDF(nc_inq_varid(group.id, "t", &group.timeId),
                                     "Time-variable lookup",
                                     "/" + group.record.name + "/t");
      if (!result)
        return result;
      if (!group.record.schedule.typeIdentifier.empty()) {
        result = detail::checkedNetCDF(
            nc_inq_varid(group.id, "portableScheduleOrdinal",
                         &group.scheduleOrdinalId),
            "Schedule-ordinal variable lookup", "/" + group.record.name);
        if (!result)
          return result;
        result = detail::checkedNetCDF(
            nc_inq_varid(group.id, "portableScheduleCursor",
                         &group.scheduleCursorId),
            "Schedule-cursor variable lookup", "/" + group.record.name);
        if (!result)
          return result;
      }
      if (group.record.containsCompleteCoefficientRestart) {
        group.coefficientReal.assign(stateLayout.coefficientFamilyCount(), -1);
        group.coefficientImag.assign(stateLayout.coefficientFamilyCount(), -1);
        for (std::size_t family = 0;
             family < stateLayout.coefficientFamilyCount(); ++family) {
          const std::string &name =
              stateLayout.coefficientFamilies()[family].identifier;
          result = detail::checkedNetCDF(
              nc_inq_varid(group.id, (name + "_real").c_str(),
                           &group.coefficientReal[family]),
              "Coefficient-variable lookup",
              "/" + group.record.name + "/" + name + "_real");
          if (!result && isDynamicsLinear) {
            result = detail::checkedNetCDF(
                nc_inq_varid(group.id, name.c_str(),
                             &group.coefficientReal[family]),
                "Initial coefficient-variable lookup",
                "/" + group.record.name + "/" + name);
            group.coefficientImag[family] = -1;
            if (!result)
              return result;
            continue;
          }
          if (!result)
            return result;
          result = detail::checkedNetCDF(
              nc_inq_varid(group.id, (name + "_imag").c_str(),
                           &group.coefficientImag[family]),
              "Coefficient-variable lookup",
              "/" + group.record.name + "/" + name + "_imag");
          if (!result)
            return result;
        }
      }
      for (auto &axisDefinition : group.unlimitedAxes) {
        result = detail::checkedNetCDF(
            nc_inq_dimid(group.id, axisDefinition.name.c_str(),
                         &axisDefinition.dimensionId),
            "Observation-axis lookup",
            "/" + group.record.name + "/" + axisDefinition.name);
        if (!result)
          return result;
        const std::string progressName =
            "portableCommitted_" + axisDefinition.name;
        result = detail::checkedNetCDF(
            nc_inq_varid(group.id, progressName.c_str(),
                         &axisDefinition.progressVariableId),
            "Observation-axis progress lookup",
            "/" + group.record.name + "/" + progressName);
        if (!result)
          return result;
      }
      for (auto &variable : group.derivedVariables) {
        const auto observerSchema = std::find_if(
            group.observerSchemas.begin(), group.observerSchemas.end(),
            [&](const auto &candidate) {
              return candidate.observerIdentifier ==
                     variable.observerIdentifier;
            });
        const bool preservesLegacyEncoding =
            observerSchema != group.observerSchemas.end() &&
            observerSchema->schema.preservesLegacyEncoding;
        if (variable.specification.scalarType ==
            WVObservationScalarType::complex64) {
          result = detail::checkedNetCDF(
              nc_inq_varid(group.id,
                           (variable.specification.name + "_real").c_str(),
                           &variable.realId),
              "Observer-variable lookup", "/" + group.record.name);
          if (!result)
            return result;
          result = detail::checkedNetCDF(
              nc_inq_varid(group.id,
                           (variable.specification.name + "_imag").c_str(),
                           &variable.imagId),
              "Observer-variable lookup", "/" + group.record.name);
        } else {
          result = detail::checkedNetCDF(
              nc_inq_varid(group.id, variable.specification.name.c_str(),
                           &variable.realId),
              "Observer-variable lookup", "/" + group.record.name);
        }
        if (!result)
          return result;
        const auto validateComponent = [&](int variableId,
                                           const std::string &variableName) {
          nc_type observedType = NC_NAT;
          auto contract = detail::checkedNetCDF(
              nc_inq_vartype(group.id, variableId, &observedType),
              "Observer-variable type inspection",
              "/" + group.record.name + "/" + variableName);
          if (!contract)
            return contract;
          const nc_type expectedType =
              variable.specification.scalarType ==
                          WVObservationScalarType::integer64
                      ? NC_INT64
                      : variable.specification.scalarType ==
                                WVObservationScalarType::boolean8
                            ? NC_UBYTE
                            : variable.specification.scalarType ==
                                      WVObservationScalarType::text
                                  ? NC_STRING
                                  : NC_DOUBLE;
          if (observedType != expectedType)
            return failure(WVCheckpointStatusCode::appendConflict,
                           "Observer-variable scalar type changed.",
                           "/" + group.record.name + "/" + variableName);
          int rank = 0;
          contract = detail::checkedNetCDF(
              nc_inq_varndims(group.id, variableId, &rank),
              "Observer-variable rank inspection",
              "/" + group.record.name + "/" + variableName);
          if (!contract)
            return contract;
          const std::size_t expectedRank =
              variable.axes.size() +
              (variable.specification.layout ==
                       WVObservationValueLayout::record
                   ? 1
                   : 0);
          if (static_cast<std::size_t>(rank) != expectedRank)
            return failure(WVCheckpointStatusCode::appendConflict,
                           "Observer-variable cadence or rank changed.",
                           "/" + group.record.name + "/" + variableName);
          std::vector<int> dimensions(expectedRank);
          if (rank > 0) {
            contract = detail::checkedNetCDF(
                nc_inq_vardimid(group.id, variableId, dimensions.data()),
                "Observer-variable dimension inspection",
                "/" + group.record.name + "/" + variableName);
            if (!contract)
              return contract;
          }
          std::vector<std::string> expectedNames;
          std::vector<const WVObservationAxis *> expectedAxes;
          if (variable.specification.layout ==
              WVObservationValueLayout::record) {
            expectedNames.push_back("t");
            expectedAxes.push_back(nullptr);
          }
          for (std::size_t reverse = variable.axes.size(); reverse > 0;
               --reverse) {
            expectedNames.push_back(variable.axes[reverse - 1].name);
            expectedAxes.push_back(&variable.axes[reverse - 1]);
          }
          for (std::size_t index = 0; index < dimensions.size(); ++index) {
            char name[NC_MAX_NAME + 1] = {};
            std::size_t length = 0;
            contract = detail::checkedNetCDF(
                nc_inq_dim(group.id, dimensions[index], name, &length),
                "Observer-variable dimension inspection",
                "/" + group.record.name + "/" + variableName);
            if (!contract)
              return contract;
            if (expectedNames[index] != name ||
                (expectedAxes[index] != nullptr &&
                 expectedAxes[index]->kind == WVObservationAxisKind::fixed &&
                 expectedAxes[index]->extent != length))
              return failure(WVCheckpointStatusCode::appendConflict,
                             "Observer-variable dimensions changed.",
                             "/" + group.record.name + "/" + variableName);
          }
          const auto requireAttribute = [&](const std::string &name,
                                            const std::string &expected,
                                            bool optional = false) {
            if (expected.empty())
              return WVCheckpointStatus::ok();
            nc_type type = NC_NAT;
            std::size_t length = 0;
            const int inquiry =
                nc_inq_att(group.id, variableId, name.c_str(), &type, &length);
            if (inquiry == NC_ENOTATT && optional)
              return WVCheckpointStatus::ok();
            auto attributeStatus = detail::checkedNetCDF(
                inquiry, "Observer-variable attribute inspection",
                "/" + group.record.name + "/" + variableName + "/@" + name);
            if (!attributeStatus)
              return attributeStatus;
            std::string observed;
            if (type == NC_CHAR) {
              observed.assign(length, '\0');
              attributeStatus = detail::checkedNetCDF(
                  nc_get_att_text(group.id, variableId, name.c_str(),
                                  observed.data()),
                  "Observer-variable attribute read",
                  "/" + group.record.name + "/" + variableName + "/@" + name);
              if (!attributeStatus)
                return attributeStatus;
            } else if (length == 1 && (expected == "0" || expected == "1")) {
              double numeric = 0.0;
              attributeStatus = detail::checkedNetCDF(
                  nc_get_att_double(group.id, variableId, name.c_str(),
                                    &numeric),
                  "Observer-variable logical-attribute read",
                  "/" + group.record.name + "/" + variableName + "/@" + name);
              if (!attributeStatus)
                return attributeStatus;
              observed = numeric == 0.0 ? "0" : numeric == 1.0 ? "1" : "";
            } else {
              return failure(WVCheckpointStatusCode::appendConflict,
                             "Observer-variable metadata type changed.",
                             "/" + group.record.name + "/" + variableName +
                                 "/@" + name);
            }
            if (observed != expected &&
                !(preservesLegacyEncoding &&
                  detail::legacyObservationAttributeMatches(
                      name, expected, observed)))
              return failure(
                  WVCheckpointStatusCode::appendConflict,
                  "Observer-variable metadata changed: expected '" + expected +
                      "' but observed '" + observed + "'.",
                  "/" + group.record.name + "/" + variableName + "/@" + name);
            return WVCheckpointStatus::ok();
          };
          contract = requireAttribute("units", variable.specification.units);
          if (contract)
            contract = requireAttribute("long_name",
                                        variable.specification.description);
          for (const auto &attribute : variable.specification.attributes) {
            if (!contract)
              break;
            contract = requireAttribute(attribute.name, attribute.value,
                                        attribute.name != "isParticle" &&
                                            attribute.name != "isTracer");
          }
          if (!contract)
            return contract;
          if ((variable.specification.layout ==
                   WVObservationValueLayout::initialValue ||
               variable.specification.layout ==
                   WVObservationValueLayout::staticValue) &&
              (variable.specification.scalarType ==
                   WVObservationScalarType::real64 ||
               variable.specification.scalarType ==
                   WVObservationScalarType::complex64)) {
            std::size_t count = 1;
            for (const auto &axis : variable.axes)
              count *= axis.extent;
            std::vector<double> values(count);
            contract = detail::checkedNetCDF(
                nc_get_var_double(group.id, variableId, values.data()),
                "Initial observer-variable read",
                "/" + group.record.name + "/" + variableName);
            if (!contract)
              return contract;
            if (std::find(values.begin(), values.end(), NC_FILL_DOUBLE) !=
                values.end())
              return failure(WVCheckpointStatusCode::incompleteRecord,
                             "Initial observer variable is unwritten.",
                             "/" + group.record.name + "/" + variableName);
          }
          return WVCheckpointStatus::ok();
        };
        result = validateComponent(variable.realId,
                                   variable.specification.scalarType ==
                                           WVObservationScalarType::complex64
                                       ? variable.specification.name + "_real"
                                       : variable.specification.name);
        if (result &&
            variable.specification.scalarType ==
                WVObservationScalarType::complex64)
          result = validateComponent(variable.imagId,
                                     variable.specification.name + "_imag");
        if (!result)
          return result;
      }
      int timeDimension = -1;
      std::size_t length = 0;
      result = detail::checkedNetCDF(
          nc_inq_dimid(group.id, "t", &timeDimension), "Time-dimension lookup",
          "/" + group.record.name + "/t");
      if (!result)
        return result;
      result = detail::checkedNetCDF(
          nc_inq_dimlen(group.id, timeDimension, &length),
          "Time-dimension length", "/" + group.record.name + "/t");
      if (!result)
        return result;
      if (length > 0) {
        std::vector<double> times(length);
        result = detail::checkedNetCDF(
            nc_get_var_double(group.id, group.timeId, times.data()),
            "Output-time read", "/" + group.record.name + "/t");
        if (!result)
          return result;
        const auto firstUncommitted =
            std::find(times.begin(), times.end(), NC_FILL_DOUBLE);
        if (firstUncommitted != times.end() &&
            std::any_of(firstUncommitted, times.end(), [](double value) {
              return value != NC_FILL_DOUBLE;
            }))
          return failure(WVCheckpointStatusCode::incompleteRecord,
                         "Output time commit markers contain an internal gap.",
                         "/" + group.record.name + "/t");
        length = static_cast<std::size_t>(
            std::distance(times.begin(), firstUncommitted));
        times.resize(length);
        group.recordCount = length;
        for (auto &axisDefinition : group.unlimitedAxes) {
          if (length == 0)
            continue;
          const std::size_t position[] = {length - 1};
          long long committedCount = -1;
          result = detail::checkedNetCDF(
              nc_get_var1_longlong(group.id,
                                   axisDefinition.progressVariableId,
                                   position, &committedCount),
              "Observation-axis progress read",
              "/" + group.record.name + "/portableCommitted_" +
                  axisDefinition.name);
          if (!result)
            return result;
          std::size_t physicalAxisLength = 0;
          result = detail::checkedNetCDF(
              nc_inq_dimlen(group.id, axisDefinition.dimensionId,
                            &physicalAxisLength),
              "Observation-axis length inspection",
              "/" + group.record.name + "/" + axisDefinition.name);
          if (!result)
            return result;
          if (committedCount < 0 ||
              static_cast<std::uint64_t>(committedCount) > physicalAxisLength)
            return failure(WVCheckpointStatusCode::incompleteRecord,
                           "Committed observation-axis count is invalid.",
                           "/" + group.record.name + "/" +
                               axisDefinition.name);
          axisDefinition.committedCount =
              static_cast<std::size_t>(committedCount);
          axisDefinition.physicalCount = physicalAxisLength;
        }
        if (group.record.schedule.typeIdentifier.empty()) {
          for (std::size_t timeIndex = 0; timeIndex < times.size();
               ++timeIndex) {
            const double observed = times[timeIndex];
            const double raw = (observed - group.record.schedule.initialTime) /
                               group.record.schedule.outputInterval;
            const auto ordinal =
                static_cast<WVOutputScheduleOrdinal>(std::llround(raw));
            const double expected = group.record.schedule.initialTime +
                                    static_cast<double>(ordinal) *
                                        group.record.schedule.outputInterval;
            const double tolerance =
                8 * std::numeric_limits<double>::epsilon() *
                std::max({1.0, std::abs(observed), std::abs(expected)});
            if (!std::isfinite(observed) ||
                std::abs(observed - expected) > tolerance || ordinal < 0 ||
                (timeIndex > 0 && !(observed > times[timeIndex - 1])))
              return failure(
                  WVCheckpointStatusCode::appendConflict,
                  "Existing output times are not a strictly increasing subset "
                  "of the configured schedule lattice.",
                  "/" + group.record.name + "/t");
            group.committedOrdinal = ordinal;
          }
        } else {
          std::vector<long long> ordinals(length);
          result = detail::checkedNetCDF(
              nc_get_var_longlong(group.id, group.scheduleOrdinalId,
                                  ordinals.data()),
              "Schedule-ordinal read", "/" + group.record.name);
          if (!result)
            return result;
          for (std::size_t timeIndex = 0; timeIndex < times.size();
               ++timeIndex) {
            if (!std::isfinite(times[timeIndex]) || ordinals[timeIndex] < 0 ||
                (timeIndex > 0 &&
                 (!(times[timeIndex] > times[timeIndex - 1]) ||
                  ordinals[timeIndex] <= ordinals[timeIndex - 1])))
              return failure(WVCheckpointStatusCode::appendConflict,
                             "Existing algorithmic output history is not "
                             "strictly monotone.",
                             "/" + group.record.name);
            const std::size_t position[] = {timeIndex};
            char *cursorText = nullptr;
            result = detail::checkedNetCDF(
                nc_get_var1_string(group.id, group.scheduleCursorId, position,
                                   &cursorText),
                "Schedule-cursor read", "/" + group.record.name);
            if (!result)
              return result;
            const std::string encoded = cursorText == nullptr ? "" : cursorText;
            if (cursorText != nullptr)
              nc_free_string(1, &cursorText);
            std::vector<std::uint8_t> cursorBytes;
            if (!hexDecode(encoded, cursorBytes) ||
                cursorBytes.size() > WVMaximumOutputScheduleCursorBytes)
              return failure(WVCheckpointStatusCode::appendConflict,
                             "Existing algorithmic output cursor is malformed.",
                             "/" + group.record.name);
            group.scheduleCursor = {};
            if (!cursorBytes.empty()) {
              const auto decoded = decodePortableTypedRecord(
                  cursorBytes, group.scheduleCursor,
                  {WVMaximumOutputScheduleCursorBytes, true, false});
              if (!decoded)
                return failure(WVCheckpointStatusCode::appendConflict,
                               decoded.message, "/" + group.record.name);
            }
            group.committedOrdinal =
                static_cast<WVOutputScheduleOrdinal>(ordinals[timeIndex]);
          }
        }

        group.hasCommittedTime = true;
        group.lastCommittedTime = times.back();

        int variableCount = 0;
        result = detail::checkedNetCDF(
            nc_inq_varids(group.id, &variableCount, nullptr),
            "Output-variable enumeration", "/" + group.record.name);
        if (!result)
          return result;
        std::vector<int> variableIds(static_cast<std::size_t>(variableCount));
        result = detail::checkedNetCDF(
            nc_inq_varids(group.id, &variableCount, variableIds.data()),
            "Output-variable enumeration", "/" + group.record.name);
        if (!result)
          return result;
        for (const int variable : variableIds) {
          nc_type type = NC_NAT;
          int dimensionCount = 0;
          result = detail::checkedNetCDF(
              nc_inq_vartype(group.id, variable, &type),
              "Output-variable type inspection", "/" + group.record.name);
          if (!result)
            return result;
          result = detail::checkedNetCDF(
              nc_inq_varndims(group.id, variable, &dimensionCount),
              "Output-variable rank inspection", "/" + group.record.name);
          if (!result)
            return result;
          if (type != NC_DOUBLE || dimensionCount == 0)
            continue;
          std::vector<int> dimensions(static_cast<std::size_t>(dimensionCount));
          result = detail::checkedNetCDF(
              nc_inq_vardimid(group.id, variable, dimensions.data()),
              "Output-variable dimension inspection", "/" + group.record.name);
          if (!result)
            return result;
          char firstName[NC_MAX_NAME + 1] = {};
          result = detail::checkedNetCDF(
              nc_inq_dimname(group.id, dimensions.front(), firstName),
              "Output-variable dimension-name inspection",
              "/" + group.record.name);
          if (!result)
            return result;
          if (std::string(firstName) != "t")
            continue;
          std::vector<std::size_t> start(dimensions.size(), 0);
          std::vector<std::size_t> count(dimensions.size(), 1);
          std::size_t slabSize = 1;
          for (std::size_t index = 1; index < dimensions.size(); ++index) {
            result = detail::checkedNetCDF(
                nc_inq_dimlen(group.id, dimensions[index], &count[index]),
                "Output-variable dimension-length inspection",
                "/" + group.record.name);
            if (!result)
              return result;
            slabSize *= count[index];
          }
          std::vector<double> slab(slabSize);
          for (std::size_t recordIndex = 0; recordIndex < length;
               ++recordIndex) {
            start.front() = recordIndex;
            result = detail::checkedNetCDF(
                nc_get_vara_double(group.id, variable, start.data(),
                                   count.data(), slab.data()),
                "Committed output-record read", "/" + group.record.name);
            if (!result)
              return result;
            if (std::find(slab.begin(), slab.end(), NC_FILL_DOUBLE) !=
                slab.end()) {
              char variableName[NC_MAX_NAME + 1] = {};
              nc_inq_varname(group.id, variable, variableName);
              return failure(WVCheckpointStatusCode::incompleteRecord,
                             "A committed output record contains unwritten "
                             "payload values.",
                             "/" + group.record.name + "/" + variableName);
            }
          }
        }
      } else {
        group.recordCount = 0;
        for (auto &axisDefinition : group.unlimitedAxes) {
          std::size_t physicalAxisLength = 0;
          result = detail::checkedNetCDF(
              nc_inq_dimlen(group.id, axisDefinition.dimensionId,
                            &physicalAxisLength),
              "Observation-axis length inspection",
              "/" + group.record.name + "/" + axisDefinition.name);
          if (!result)
            return result;
          axisDefinition.physicalCount = physicalAxisLength;
          if (physicalAxisLength != 0)
            return failure(WVCheckpointStatusCode::incompleteRecord,
                           "An empty output group has committed ragged storage.",
                           "/" + group.record.name + "/" +
                               axisDefinition.name);
        }
      }
    }
    return WVCheckpointStatus::ok();
  }

  void rebuildProgress() {
    destinationProgress.clear();
    for (const auto &file : files)
      for (const auto &group : file.groups) {
        WVOutputDestinationProgress progress;
        progress.fileIdentifier = file.record.identifier;
        progress.groupIdentifier = group.record.identifier;
        progress.recordCount = group.recordCount;
        progress.hasCommittedTime = group.hasCommittedTime;
        progress.lastCommittedTime = group.lastCommittedTime;
        progress.committedScheduleCursor = {group.committedOrdinal,
                                            group.scheduleCursor};
        progress.unlimitedAxes.reserve(group.unlimitedAxes.size());
        for (const auto &axis : group.unlimitedAxes)
          progress.unlimitedAxes.push_back(
              {axis.name, axis.committedCount, axis.physicalCount});
        destinationProgress.push_back(std::move(progress));
      }
  }

  void updateRetainedStorageMetric() {
    std::size_t bytes = sizeof(*this) + stateLayout.persistentBytes() +
                        (transformAdapter == nullptr
                             ? 0
                             : transformAdapter->persistentBytes()) +
                        files.capacity() * sizeof(File) +
                        destinationProgress.capacity() *
                            sizeof(WVOutputDestinationProgress) +
                        committedPreparedRoutes.capacity() *
                            sizeof(std::pair<std::size_t, std::size_t>) +
                        preparedObservationStatus.message.capacity() +
                        preparedObservationStatus.location.capacity();
    for (const auto &file : files) {
      bytes += file.destination.native().capacity() *
                   sizeof(std::filesystem::path::value_type) +
               file.staging.native().capacity() *
                   sizeof(std::filesystem::path::value_type) +
               file.groups.capacity() * sizeof(Group);
      for (const auto &group : file.groups) {
        bytes += group.scheduleCursor.persistentBytes() -
                 sizeof(WVPortableTypedRecord) +
                 group.coefficientReal.capacity() * sizeof(int) +
                 group.coefficientImag.capacity() * sizeof(int) +
                 group.derivedVariables.capacity() * sizeof(Variable) +
                 group.observerSchemas.capacity() * sizeof(ObserverSchema) +
                 group.unlimitedAxes.capacity() * sizeof(UnlimitedAxis) +
                 group.preparedBatchIndices.capacity() * sizeof(std::size_t) +
                 group.preparedUnlimitedAxisIncrements.capacity() *
                     sizeof(std::size_t) +
                 group.realWriteScratch.capacity() * sizeof(double) +
                 group.imaginaryWriteScratch.capacity() * sizeof(double) +
                 group.integerWriteScratch.capacity() * sizeof(long long) +
                 group.textWriteScratch.capacity() * sizeof(const char *);
        for (const auto &variable : group.derivedVariables)
          bytes += variable.observerIdentifier.capacity() +
                   variable.schemaIdentifier.capacity() +
                   variable.specification.identifier.capacity() +
                   variable.specification.name.capacity() +
                   variable.specification.units.capacity() +
                   variable.specification.description.capacity() +
                   variable.specification.raggedChildAxisIdentifier.capacity() +
                   variable.specification.dimensionIdentifiers.capacity() *
                       sizeof(std::string) +
                   variable.axes.capacity() * sizeof(WVObservationAxis) +
                   variable.unlimitedAxisIndices.capacity() *
                       sizeof(std::size_t) +
                   variable.writeStart.capacity() * sizeof(std::size_t) +
                   variable.writeCount.capacity() * sizeof(std::size_t) +
                   variable.specification.attributes.capacity() *
                       sizeof(WVObservationAttribute);
        for (const auto &variable : group.derivedVariables) {
          for (const auto &name : variable.specification.dimensionIdentifiers)
            bytes += name.capacity();
          for (const auto &axisDefinition : variable.axes)
            bytes += axisDefinition.identifier.capacity() +
                     axisDefinition.name.capacity();
          for (const auto &attribute : variable.specification.attributes)
            bytes += attribute.name.capacity() + attribute.value.capacity();
        }
        for (const auto &axisDefinition : group.unlimitedAxes)
          bytes += axisDefinition.name.capacity();
        for (const auto &observerSchema : group.observerSchemas) {
          bytes += observerSchema.observerIdentifier.capacity() +
                   observationSchemaRetainedBytes(observerSchema.schema) +
                   observerSchema.variableAxisIndices.capacity() *
                       sizeof(std::vector<std::size_t>) +
                   observerSchema.persistenceVariableIndices.capacity() *
                       sizeof(std::size_t) +
                   observerSchema.unlimitedAxisIndices.capacity() *
                       sizeof(std::size_t) +
                   observerSchema.raggedChildAxisIndices.capacity() *
                       sizeof(std::size_t) +
                   observerSchema.preparedValueIndices.capacity() *
                       sizeof(std::size_t) +
                   observerSchema.preparedUnlimitedExtents.capacity() *
                       sizeof(std::size_t) +
                   observerSchema.preparedObservedVariables.capacity() *
                       sizeof(std::uint8_t);
          for (const auto &axisIndices :
               observerSchema.variableAxisIndices)
            bytes += axisIndices.capacity() * sizeof(std::size_t);
        }
      }
    }
    for (const auto &item : destinationProgress) {
      bytes +=
          item.fileIdentifier.capacity() + item.groupIdentifier.capacity() +
          item.committedScheduleCursor.values.persistentBytes() -
              sizeof(WVPortableTypedRecord) +
          item.unlimitedAxes.capacity() *
              sizeof(WVOutputDestinationAxisProgress);
      for (const auto &axis : item.unlimitedAxes)
        bytes += axis.axisIdentifier.capacity();
    }
    bytes += preparedObservationBatches.capacity() *
             sizeof(PreparedObservationBatch);
    for (const auto &prepared : preparedObservationBatches)
      bytes += prepared.batch.metrics().retainedStorageBytes;
    metrics.retainedStorageBytes = bytes;
  }

  WVCheckpointStatus openExistingFiles(
      const std::vector<WVOutputDestinationProgress> &expectedProgress) {
    files.clear();
    files.reserve(descriptor.record().outputFiles.size());
    for (const auto &record : descriptor.record().outputFiles) {
      WVCheckpointInspection inspection;
      auto result = WVCheckpointReader::inspect(record.destination,
                                                *catalog,
                                                inspection);
      if (!result)
        return result;
      if (transformAdapter == nullptr ||
          !transformAdapter->sameConfiguration(inspection))
        return failure(WVCheckpointStatusCode::schemaMismatch,
                       "Append model configuration does not match the file.",
                       record.destination);
      File file(record);
      file.destination =
          std::filesystem::absolute(record.destination).lexically_normal();
      int id = -1;
      result = detail::checkedNetCDF(
          nc_open(file.destination.c_str(), NC_NOWRITE, &id),
          "Append output-file read-only preflight open",
          file.destination.string());
      if (!result)
        return result;
      file.id = id;
      files.push_back(std::move(file));
      auto &openFile = files.back();
      openFile.groups.reserve(record.groups.size());
      for (const auto &groupRecord : record.groups) {
        Group group(groupRecord);
        // Rebuild variable contracts without mutating the existing file.
        for (const auto &identifier : groupRecord.observerIdentifiers) {
          const auto *recordPointer = observer(identifier);
          if (recordPointer == nullptr)
            return failure(WVCheckpointStatusCode::schemaMismatch,
                           "Append group references an unknown observer.",
                           openFile.destination.string());
          if (sampleSource != nullptr) {
            WVObservationSchema schema;
            const auto sourceStatus = sampleSource->observationSchema(
                *recordPointer, schema);
            if (!sourceStatus)
              return failure(WVCheckpointStatusCode::unsupportedObserver,
                             sourceStatus.message,
                             openFile.destination.string());
            const auto schemaStatus = validateObservationSchema(schema);
            if (!schemaStatus)
              return failure(WVCheckpointStatusCode::schemaMismatch,
                             schemaStatus.message,
                             openFile.destination.string());
            for (const auto &specification : schema.variables) {
              const bool canonicalCoefficient =
                  groupRecord.containsCompleteCoefficientRestart &&
                  coefficientMetadata(specification.name) != nullptr;
              if (!canonicalCoefficient) {
                const auto existing = std::find_if(
                    group.derivedVariables.begin(),
                    group.derivedVariables.end(), [&](const auto &variable) {
                      return variable.specification.name == specification.name;
                    });
                if (existing == group.derivedVariables.end()) {
                  Variable variable;
                  variable.observerIdentifier = identifier;
                  variable.schemaIdentifier = schema.identifier;
                  variable.schemaVersion = schema.version;
                  variable.specification = specification;
                  for (const auto &axisIdentifier :
                       specification.dimensionIdentifiers) {
                    const auto *axisDefinition = axis(schema, axisIdentifier);
                    if (axisDefinition == nullptr)
                      return failure(
                          WVCheckpointStatusCode::schemaMismatch,
                          "Observation variable references an unknown axis.",
                          openFile.destination.string());
                    variable.axes.push_back(*axisDefinition);
                    if (axisDefinition->kind ==
                            WVObservationAxisKind::unlimited &&
                        std::none_of(
                            group.unlimitedAxes.begin(),
                            group.unlimitedAxes.end(), [&](const auto &item) {
                              return item.name == axisDefinition->name;
                            })) {
                      UnlimitedAxis unlimited;
                      unlimited.name = axisDefinition->name;
                      group.unlimitedAxes.push_back(std::move(unlimited));
                    }
                  }
                  group.derivedVariables.push_back(std::move(variable));
                } else
                  return failure(
                      WVCheckpointStatusCode::appendConflict,
                      "Observers declare duplicate persisted variable names "
                      "without an alias contract.",
                      openFile.destination.string());
              }
            }
            group.observerSchemas.push_back(
                {identifier, std::move(schema)});
          }
        }
        openFile.groups.push_back(std::move(group));
      }
      result = resolveFile(openFile);
      if (!result)
        return result;
    }

    rebuildProgress();
    if (!sameDestinationProgress(destinationProgress, expectedProgress))
      return failure(WVCheckpointStatusCode::appendConflict,
                     "Append destination progress changed after compilation.",
                     "/");

    // Only after the entire destination set has passed read-only graph,
    // schema, cursor, record, time-last, payload, and ragged-offset validation
    // do any files reopen with mutation capability.
    for (auto &file : files) {
      const int id = std::exchange(file.id, -1);
      const int code = nc_close(id);
      if (code != NC_NOERR)
        return detail::netcdfFailure(code, "Append preflight close",
                                     file.destination.string());
    }
    for (auto &file : files) {
      int id = -1;
      auto result = detail::checkedNetCDF(
          nc_open(file.destination.c_str(), NC_WRITE, &id),
          "Append output-file mutation open", file.destination.string());
      if (!result)
        return result;
      file.id = id;
      result = resolveFile(file);
      if (!result)
        return result;
    }
    return WVCheckpointStatus::ok();
  }

  WVCheckpointStatus writeRealSlab(int group, int variable,
                                   std::size_t recordIndex,
                                   const double *values,
                                   std::size_t elementCount,
                                   const std::string &path) {
    int dimensionCount = 0;
    auto result =
        detail::checkedNetCDF(nc_inq_varndims(group, variable, &dimensionCount),
                              "Output-variable rank inspection", path);
    if (!result)
      return result;
    std::vector<int> dimensionIds(static_cast<std::size_t>(dimensionCount));
    result = detail::checkedNetCDF(
        nc_inq_vardimid(group, variable, dimensionIds.data()),
        "Output-variable dimension inspection", path);
    if (!result)
      return result;
    std::vector<std::size_t> start(static_cast<std::size_t>(dimensionCount), 0);
    std::vector<std::size_t> count(static_cast<std::size_t>(dimensionCount), 1);
    start.front() = recordIndex;
    std::size_t observed = 1;
    for (std::size_t index = 1; index < dimensionIds.size(); ++index) {
      result = detail::checkedNetCDF(
          nc_inq_dimlen(group, dimensionIds[index], &count[index]),
          "Output-variable dimension length", path);
      if (!result)
        return result;
      observed *= count[index];
    }
    if (observed != elementCount)
      return failure(WVCheckpointStatusCode::shapeMismatch,
                     "Output value and NetCDF slab sizes differ.", path);
    return detail::checkedNetCDF(
        nc_put_vara_double(group, variable, start.data(), count.data(), values),
        "Output-variable write", path);
  }

  WVCheckpointStatus writeComplexSlab(int group, int realId, int imagId,
                                      std::size_t recordIndex,
                                      const WVComplex64 *values,
                                      std::size_t elementCount,
                                      const std::string &path) {
    std::vector<double> real(elementCount);
    std::vector<double> imag(elementCount);
    for (std::size_t index = 0; index < elementCount; ++index) {
      real[index] = values[index].real;
      imag[index] = values[index].imag;
    }
    auto result = writeRealSlab(group, realId, recordIndex, real.data(),
                                elementCount, path + "_real");
    if (!result)
      return result;
    return writeRealSlab(group, imagId, recordIndex, imag.data(), elementCount,
                         path + "_imag");
  }

  WVCheckpointStatus writeRoute(const WVOutputEvent &event,
                                const WVOutputRouteView &route,
                                WVOutputDeliveryResult &delivery) {
    const auto payloadStarted = std::chrono::steady_clock::now();
    if (route.fileOrdinal >= files.size() ||
        route.groupOrdinal >= files[route.fileOrdinal].groups.size())
      return failure(WVCheckpointStatusCode::schemaMismatch,
                     "Output route is outside the configured file graph.", "/");
    auto &file = files[route.fileOrdinal];
    auto &group = file.groups[route.groupOrdinal];
    const std::size_t index = group.recordCount;
    constexpr auto unresolved = WVNoResolvedObservationVariable;
    std::fill(group.preparedUnlimitedAxisIncrements.begin(),
              group.preparedUnlimitedAxisIncrements.end(), unresolved);
    if (group.preparedBatchIndices.size() != group.observerSchemas.size())
      return failure(WVCheckpointStatusCode::unsupportedObserver,
                     "The route has no compiled occurrence-batch bindings.",
                     "/" + group.record.name);
    for (std::size_t observerIndex = 0;
         observerIndex < group.observerSchemas.size(); ++observerIndex) {
      const auto &observerSchema = group.observerSchemas[observerIndex];
      const auto batchIndex = group.preparedBatchIndices[observerIndex];
      const auto *batch = batchIndex < preparedObservationBatches.size()
                              ? &preparedObservationBatches[batchIndex].batch
                              : nullptr;
      if (batch == nullptr)
        return failure(WVCheckpointStatusCode::unsupportedObserver,
                       "The exact event has no staged observation batch.",
                       "/" + group.record.name);
      for (const auto &value : batch->values) {
        const auto variableIndex = value.resolvedVariableIndex;
        if (variableIndex >= observerSchema.schema.variables.size())
          return failure(WVCheckpointStatusCode::schemaMismatch,
                         "Observation batch has no resolved value slot.",
                         "/" + group.record.name);
        for (std::size_t dimension = 0;
             dimension <
             observerSchema.variableAxisIndices[variableIndex].size();
             ++dimension) {
          const auto schemaAxisIndex =
              observerSchema.variableAxisIndices[variableIndex][dimension];
          const auto &axisDefinition =
              observerSchema.schema.axes[schemaAxisIndex];
          if (axisDefinition.kind != WVObservationAxisKind::unlimited)
            continue;
          const auto groupAxisIndex =
              observerSchema.unlimitedAxisIndices[schemaAxisIndex];
          if (groupAxisIndex >=
              group.preparedUnlimitedAxisIncrements.size())
            return failure(
                WVCheckpointStatusCode::schemaMismatch,
                "Observation batch has an unresolved unlimited axis.",
                "/" + group.record.name);
          auto &increment =
              group.preparedUnlimitedAxisIncrements[groupAxisIndex];
          if (increment == unresolved)
            increment = value.extents[dimension];
          else if (increment != value.extents[dimension])
            return failure(
                WVCheckpointStatusCode::shapeMismatch,
                "Observation batches disagree on a shared unlimited axis.",
                "/" + group.record.name + "/" + axisDefinition.name);
        }
      }
    }
    if (group.record.containsCompleteCoefficientRestart &&
        !isDynamicsLinear) {
      for (std::size_t family = 0;
           family < stateLayout.coefficientFamilyCount(); ++family) {
        const auto &familyLayout = stateLayout.coefficientFamilies()[family];
        const auto values = coefficientFamilyView(stateLayout, event.state,
                                                  family);
        if (values.data == nullptr)
          return failure(WVCheckpointStatusCode::shapeMismatch,
                         "The event omits a transform coefficient family.",
                         "/" + group.record.name + "/" +
                             familyLayout.identifier);
        auto result = writeComplexSlab(
            group.id, group.coefficientReal[family],
            group.coefficientImag[family], index, values.data,
            familyLayout.elementCount,
            "/" + group.record.name + "/" + familyLayout.identifier);
        if (!result)
          return result;
        delivery.writeCount += 2;
        delivery.writtenBytes +=
            2 * familyLayout.elementCount * sizeof(double);
      }
    }
    for (std::size_t schemaIndex = 0;
         schemaIndex < group.observerSchemas.size(); ++schemaIndex) {
      const auto &observerSchema = group.observerSchemas[schemaIndex];
      const auto batchIndex = group.preparedBatchIndices[schemaIndex];
      if (batchIndex >= preparedObservationBatches.size())
        return failure(WVCheckpointStatusCode::unsupportedObserver,
                       "The occurrence-batch binding is invalid.",
                       "/" + group.record.name);
      const auto &batch = preparedObservationBatches[batchIndex].batch;
      for (const auto &value : batch.values) {
        const auto variableIndex = value.resolvedVariableIndex;
        if (variableIndex >= observerSchema.persistenceVariableIndices.size())
          return failure(WVCheckpointStatusCode::schemaMismatch,
                         "Observation value has no resolved persistence slot.",
                         "/" + group.record.name);
        const auto persistenceIndex =
            observerSchema.persistenceVariableIndices[variableIndex];
        if (persistenceIndex == unresolved)
          continue;
        if (persistenceIndex >= group.derivedVariables.size())
          return failure(WVCheckpointStatusCode::schemaMismatch,
                         "Observation persistence binding is invalid.",
                         "/" + group.record.name);
        auto &variable = group.derivedVariables[persistenceIndex];
        auto result = writeObservationValue(
            group, variable, value, index,
            "/" + group.record.name + "/" + variable.specification.name);
        if (!result)
          return result;
        if (value.elementCount() > 0)
          delivery.writeCount +=
              value.scalarType == WVObservationScalarType::complex64 ? 2 : 1;
        delivery.writtenBytes += value.liveBytes();
      }
    }
    const std::size_t start[] = {index};
    const std::size_t count[] = {1};
    for (std::size_t axisIndex = 0; axisIndex < group.unlimitedAxes.size();
         ++axisIndex) {
      const auto &axisDefinition = group.unlimitedAxes[axisIndex];
      const auto increment = group.preparedUnlimitedAxisIncrements[axisIndex];
      const std::size_t value =
          axisDefinition.committedCount +
          (increment == unresolved ? 0 : increment);
      if (value > static_cast<std::size_t>(
                      std::numeric_limits<long long>::max()) ||
          value < axisDefinition.committedCount)
        return failure(WVCheckpointStatusCode::invalidValue,
                       "Committed observation-axis count overflowed.",
                       "/" + group.record.name + "/" + axisDefinition.name);
      const auto persisted = static_cast<long long>(value);
      auto progressResult = detail::checkedNetCDF(
          nc_put_var1_longlong(group.id, axisDefinition.progressVariableId,
                               start, &persisted),
          "Observation-axis progress write",
          "/" + group.record.name + "/portableCommitted_" +
              axisDefinition.name);
      if (!progressResult)
        return progressResult;
      ++delivery.writeCount;
      delivery.writtenBytes += sizeof(std::int64_t);
    }
    if (!group.record.schedule.typeIdentifier.empty()) {
      if (route.proposedScheduleCursor == nullptr)
        return failure(WVCheckpointStatusCode::invalidValue,
                       "Algorithmic output route has no proposed cursor.",
                       "/" + group.record.name);
      std::vector<std::uint8_t> cursorBytes;
      if (!isCanonicalEmptyPortableTypedRecord(
              *route.proposedScheduleCursor)) {
        auto cursorStatus = encodePortableTypedRecord(
            *route.proposedScheduleCursor, cursorBytes);
        if (!cursorStatus)
          return failure(WVCheckpointStatusCode::invalidValue,
                         cursorStatus.message, "/" + group.record.name);
      }
      if (cursorBytes.size() > WVMaximumOutputScheduleCursorBytes)
        return failure(WVCheckpointStatusCode::invalidValue,
                       "Algorithmic output cursor exceeds 4 KiB.",
                       "/" + group.record.name);
      const auto cursorHex = hexEncode(cursorBytes);
      const char *cursorText = cursorHex.c_str();
      const auto ordinal = static_cast<long long>(route.scheduleOrdinal);
      auto scheduleResult = detail::checkedNetCDF(
          nc_put_var1_longlong(group.id, group.scheduleOrdinalId, start,
                               &ordinal),
          "Schedule-ordinal write", "/" + group.record.name);
      if (!scheduleResult)
        return scheduleResult;
      scheduleResult = detail::checkedNetCDF(
          nc_put_var1_string(group.id, group.scheduleCursorId, start,
                             &cursorText),
          "Schedule-cursor write", "/" + group.record.name);
      if (!scheduleResult)
        return scheduleResult;
      delivery.writeCount += 2;
      delivery.writtenBytes += sizeof(std::int64_t) + cursorHex.size();
    }
    auto result = detail::checkedNetCDF(
        nc_put_vara_double(group.id, group.timeId, start, count,
                           &event.scheduledTime),
        "Output-time commit", "/" + group.record.name + "/t");
    if (!result)
      return result;
    metrics.payloadWriteSeconds +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      payloadStarted)
            .count();
    const auto synchronizationStarted = std::chrono::steady_clock::now();
    result = detail::checkedNetCDF(nc_sync(file.id), "Output-route sync",
                                   file.destination.string());
    if (!result)
      return result;
    metrics.synchronizationSeconds +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      synchronizationStarted)
            .count();
    ++group.recordCount;
    group.hasCommittedTime = true;
    group.lastCommittedTime = event.scheduledTime;
    for (std::size_t axisIndex = 0; axisIndex < group.unlimitedAxes.size();
         ++axisIndex) {
      auto &axisDefinition = group.unlimitedAxes[axisIndex];
      const auto increment = group.preparedUnlimitedAxisIncrements[axisIndex];
      if (increment != unresolved) {
        axisDefinition.committedCount += increment;
        axisDefinition.physicalCount = axisDefinition.committedCount;
      }
    }
    group.committedOrdinal = route.scheduleOrdinal;
    if (route.proposedScheduleCursor != nullptr)
      group.scheduleCursor = *route.proposedScheduleCursor;
    ++metrics.committedRecordCount;
    ++metrics.synchronizationCount;
    metrics.writtenBytes += delivery.writtenBytes + sizeof(double);
    delivery.writeCount += 1;
    delivery.writtenBytes += sizeof(double);
    rebuildProgress();
    return WVCheckpointStatus::ok();
  }

  WVCheckpointStatus close() noexcept {
    WVCheckpointStatus first = WVCheckpointStatus::ok();
    for (auto &file : files) {
      if (file.id >= 0) {
        const int id = std::exchange(file.id, -1);
        const int code = nc_close(id);
        if (code != NC_NOERR && first)
          first = detail::netcdfFailure(code, "Output-file close",
                                        file.destination.string());
      }
      if (!file.staging.empty()) {
        std::error_code ignored;
        std::filesystem::remove(file.staging, ignored);
      }
    }
    closed = true;
    return first;
  }
};

WVModelOutputNetCDFSink::WVModelOutputNetCDFSink() : impl_(new Impl) {}
WVModelOutputNetCDFSink::~WVModelOutputNetCDFSink() = default;
WVModelOutputNetCDFSink::WVModelOutputNetCDFSink(
    WVModelOutputNetCDFSink &&) noexcept = default;
WVModelOutputNetCDFSink &WVModelOutputNetCDFSink::operator=(
    WVModelOutputNetCDFSink &&) noexcept = default;

WVCheckpointStatus WVModelOutputNetCDFSink::createNew(
    const WVModelOutputNetCDFConfiguration &configuration,
    const WVPortableObserverDescriptor &descriptor, const WVOutputPlan &plan,
    const WVIntegrationStateLayout &stateLayout,
    WVObserverSampleSource *sampleSource, WVModelOutputNetCDFSink &sink) {
  if (sampleSource != nullptr) {
    const auto status = sampleSource->preflight(plan);
    if (!status)
      return failure(WVCheckpointStatusCode::unsupportedObserver,
                     status.message, "/observingSystems");
  }
  auto catalogStatus = validateCatalogIdentity(configuration, descriptor);
  if (!catalogStatus)
    return catalogStatus;
  try {
    auto candidate = std::make_unique<Impl>();
    candidate->catalog = configuration.catalog;
    auto result = detail::createModelOutputTransformAdapter(
        configuration.checkpointTemplate, candidate->transformAdapter);
    if (!result)
      return result;
    candidate->constructionCheckpoint = &configuration.checkpointTemplate;
    candidate->isDynamicsLinear = configuration.isDynamicsLinear;
    candidate->descriptor = descriptor;
    candidate->stateLayout = stateLayout;
    candidate->sampleSource = sampleSource;
    result = candidate->validateConfiguration();
    if (!result)
      return result;
    const auto planStatus = candidate->validateAndCompilePlan(plan);
    if (!planStatus)
      return failure(WVCheckpointStatusCode::schemaMismatch, planStatus.message,
                     "/outputPlan");
    result = candidate->validateDestinations(false);
    if (!result)
      return result;
    result = candidate->stageFiles();
    if (!result)
      return result;
    result = candidate->commitStagedFiles();
    if (!result)
      return result;
    candidate->rebuildProgress();
    candidate->metrics.fileCount = candidate->files.size();
    for (const auto &file : candidate->files)
      candidate->metrics.groupCount += file.groups.size();
    candidate->metrics.initializedFileCount = candidate->files.size();
    candidate->constructionCheckpoint = nullptr;
    candidate->updateRetainedStorageMetric();
    sink.impl_ = std::move(candidate);
    return WVCheckpointStatus::ok();
  } catch (const std::bad_alloc &) {
    return failure(WVCheckpointStatusCode::writeFailure,
                   "Output persistence allocation failed.", "/");
  } catch (const std::exception &exception) {
    return failure(WVCheckpointStatusCode::writeFailure,
                   "Output persistence creation failed: " +
                       std::string(exception.what()),
                   "/");
  }
}

WVCheckpointStatus WVModelOutputNetCDFSink::replaceExisting(
    const WVModelOutputNetCDFConfiguration &configuration,
    const WVPortableObserverDescriptor &descriptor, const WVOutputPlan &plan,
    const WVIntegrationStateLayout &stateLayout,
    WVObserverSampleSource *sampleSource, WVModelOutputNetCDFSink &sink) {
  if (sampleSource != nullptr) {
    const auto status = sampleSource->preflight(plan);
    if (!status)
      return failure(WVCheckpointStatusCode::unsupportedObserver,
                     status.message, "/observingSystems");
  }
  auto catalogStatus = validateCatalogIdentity(configuration, descriptor);
  if (!catalogStatus)
    return catalogStatus;
  try {
    auto candidate = std::make_unique<Impl>();
    candidate->catalog = configuration.catalog;
    auto result = detail::createModelOutputTransformAdapter(
        configuration.checkpointTemplate, candidate->transformAdapter);
    if (!result)
      return result;
    candidate->constructionCheckpoint = &configuration.checkpointTemplate;
    candidate->isDynamicsLinear = configuration.isDynamicsLinear;
    candidate->descriptor = descriptor;
    candidate->stateLayout = stateLayout;
    candidate->sampleSource = sampleSource;
    result = candidate->validateConfiguration();
    if (!result)
      return result;
    const auto planStatus = candidate->validateAndCompilePlan(plan);
    if (!planStatus)
      return failure(WVCheckpointStatusCode::schemaMismatch, planStatus.message,
                     "/outputPlan");
    result = candidate->validateDestinations(true);
    if (!result)
      return result;
    result = candidate->stageFiles();
    if (!result)
      return result;
    result = candidate->replaceStagedFiles();
    if (!result)
      return result;
    candidate->rebuildProgress();
    candidate->metrics.fileCount = candidate->files.size();
    for (const auto &file : candidate->files)
      candidate->metrics.groupCount += file.groups.size();
    candidate->metrics.initializedFileCount = candidate->files.size();
    candidate->constructionCheckpoint = nullptr;
    candidate->updateRetainedStorageMetric();
    sink.impl_ = std::move(candidate);
    return WVCheckpointStatus::ok();
  } catch (const std::bad_alloc &) {
    return failure(WVCheckpointStatusCode::writeFailure,
                   "Output replacement allocation failed.", "/");
  } catch (const std::exception &exception) {
    return failure(
        WVCheckpointStatusCode::writeFailure,
        "Output replacement failed: " + std::string(exception.what()), "/");
  }
}

WVCheckpointStatus WVModelOutputNetCDFSink::openAppend(
    const WVModelOutputNetCDFConfiguration &configuration,
    const WVPortableObserverDescriptor &descriptor, const WVOutputPlan &plan,
    const WVIntegrationStateLayout &stateLayout,
    WVObserverSampleSource *sampleSource,
    const std::vector<WVOutputDestinationProgress> &expectedDestinationProgress,
    WVModelOutputNetCDFSink &sink) {
  if (sampleSource != nullptr) {
    const auto status = sampleSource->preflight(plan);
    if (!status)
      return failure(WVCheckpointStatusCode::unsupportedObserver,
                     status.message, "/observingSystems");
  }
  auto catalogStatus = validateCatalogIdentity(configuration, descriptor);
  if (!catalogStatus)
    return catalogStatus;
  try {
    auto candidate = std::make_unique<Impl>();
    candidate->catalog = configuration.catalog;
    auto result = detail::createModelOutputTransformAdapter(
        configuration.checkpointTemplate, candidate->transformAdapter);
    if (!result)
      return result;
    candidate->constructionCheckpoint = &configuration.checkpointTemplate;
    candidate->isDynamicsLinear = configuration.isDynamicsLinear;
    candidate->descriptor = descriptor;
    candidate->stateLayout = stateLayout;
    candidate->sampleSource = sampleSource;
    candidate->appendMode = true;
    result = candidate->validateConfiguration();
    if (!result)
      return result;
    const auto planStatus = candidate->validateAndCompilePlan(plan);
    if (!planStatus)
      return failure(WVCheckpointStatusCode::schemaMismatch, planStatus.message,
                     "/outputPlan");
    result = candidate->validateDestinations(true);
    if (!result)
      return result;
    result = candidate->openExistingFiles(expectedDestinationProgress);
    if (!result)
      return result;
    candidate->metrics.fileCount = candidate->files.size();
    for (const auto &file : candidate->files)
      candidate->metrics.groupCount += file.groups.size();
    candidate->constructionCheckpoint = nullptr;
    candidate->updateRetainedStorageMetric();
    sink.impl_ = std::move(candidate);
    return WVCheckpointStatus::ok();
  } catch (const std::exception &exception) {
    return failure(WVCheckpointStatusCode::openFailure,
                   "Output append failed: " + std::string(exception.what()),
                   "/");
  }
}

WVKernelStatus WVModelOutputNetCDFSink::preflight(const WVOutputPlan &plan) {
  if (!impl_ || impl_->closed)
    return {WVKernelStatusCode::invalidConfiguration,
            "NetCDF output sink is closed."};
  if (impl_->sampleSource != nullptr) {
    const auto sourceStatus = impl_->sampleSource->preflight(plan);
    if (!sourceStatus)
      return sourceStatus;
  }
  const auto planStatus = impl_->validateAndCompilePlan(plan);
  if (!planStatus)
    return planStatus;
  impl_->updateRetainedStorageMetric();
  impl_->preflightComplete = true;
  return WVKernelStatus::ok();
}

WVKernelStatus
WVModelOutputNetCDFSink::deliver(const WVOutputEvent &event,
                                 const WVOutputRouteView &route,
                                 WVOutputDeliveryResult &result) {
  if (!impl_ || !impl_->preflightComplete || impl_->closed)
    return {WVKernelStatusCode::invalidConfiguration,
            "NetCDF output sink has not completed preflight."};
  const auto preparation = impl_->prepareObservationBatches(event);
  if (!preparation) {
    ++impl_->metrics.failureCount;
    return kernelFailure(preparation);
  }
  const auto status = impl_->writeRoute(event, route, result);
  if (!status) {
    // A failed NetCDF call can occur after a conversion scratch buffer grows.
    // Publish that still-retained capacity before returning so retry and
    // maximum-live ledgers do not silently omit it.
    impl_->updatePreparedObservationMetrics();
    ++impl_->metrics.failureCount;
    return kernelFailure(status);
  }
  impl_->completePreparedRoute(event, route);
  return WVKernelStatus::ok();
}

const std::vector<WVOutputDestinationProgress> &
WVModelOutputNetCDFSink::destinationProgress() const noexcept {
  static const std::vector<WVOutputDestinationProgress> empty;
  return impl_ ? impl_->destinationProgress : empty;
}

const WVModelOutputNetCDFMetrics &
WVModelOutputNetCDFSink::metrics() const noexcept {
  static const WVModelOutputNetCDFMetrics empty;
  return impl_ ? impl_->metrics : empty;
}

WVCheckpointStatus WVModelOutputNetCDFSink::close() noexcept {
  return impl_ ? impl_->close() : WVCheckpointStatus::ok();
}

} // namespace wavevortex::runtime
