#include "WaveVortexRuntime/WVModelOutputConfiguration.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVObserverOutputProvider.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <utility>

namespace wavevortex::runtime {
namespace {

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

WVKernelStatus fromCheckpoint(const WVCheckpointStatus &status) {
  if (status)
    return WVKernelStatus::ok();
  return {status.code == WVCheckpointStatusCode::unsupportedObserver
              ? WVKernelStatusCode::unsupportedOperation
              : WVKernelStatusCode::invalidConfiguration,
          status.message + (status.location.empty() ? "" :
                            " [" + status.location + "]")};
}

std::string normalizedDestination(const std::string &value) {
  return std::filesystem::absolute(std::filesystem::path(value))
      .lexically_normal()
      .string();
}

std::string stableIdentifier(const std::string &prefix,
                             const std::string &value) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  for (const auto character : value) {
    hash ^= static_cast<unsigned char>(character);
    hash *= UINT64_C(1099511628211);
  }
  std::ostringstream stream;
  stream << prefix << '-' << std::hex << std::setw(16) << std::setfill('0')
         << hash;
  return stream.str();
}

bool sameTypedRecord(const WVPortableTypedRecord &left,
                     const WVPortableTypedRecord &right) noexcept;

bool sameSchedule(const WVOutputScheduleRecord &left,
                  const WVOutputScheduleRecord &right) noexcept {
  if (left.outputInterval != right.outputInterval ||
      left.initialTime != right.initialTime ||
      left.finalTime != right.finalTime ||
      left.typeIdentifier != right.typeIdentifier ||
      left.contractVersion != right.contractVersion ||
      left.configuration.schemaIdentifier !=
          right.configuration.schemaIdentifier ||
      left.configuration.schemaVersion != right.configuration.schemaVersion ||
      left.configuration.values.size() != right.configuration.values.size())
    return false;
  for (std::size_t index = 0; index < left.configuration.values.size(); ++index) {
    const auto &a = left.configuration.values[index];
    const auto &b = right.configuration.values[index];
    if (a.name != b.name || a.dimensions != b.dimensions ||
        a.storage != b.storage)
      return false;
  }
  return true;
}

bool sameStateBlock(const WVStateBlockRecord &left,
                    const WVStateBlockRecord &right) noexcept {
  return left.identifier == right.identifier &&
         left.scalarType == right.scalarType &&
         left.dimensions == right.dimensions &&
         left.toleranceKind == right.toleranceKind &&
         left.absoluteTolerance == right.absoluteTolerance &&
         left.ownership == right.ownership &&
         left.restartRequirement == right.restartRequirement;
}

bool sameObserver(const WVObserverRecord &left,
                  const WVObserverRecord &right) noexcept {
  return left.identifier == right.identifier && left.name == right.name &&
         left.typeIdentifier == right.typeIdentifier &&
         left.contractVersion == right.contractVersion &&
         sameTypedRecord(left.configuration, right.configuration) &&
         left.stateBlockIdentifiers == right.stateBlockIdentifiers &&
         left.fieldNames == right.fieldNames && left.x == right.x &&
         left.y == right.y && left.z == right.z &&
         left.isXYOnly == right.isXYOnly &&
         left.shouldAntialias == right.shouldAntialias &&
         left.advectionInterpolation == right.advectionInterpolation &&
         left.trackedFieldInterpolation == right.trackedFieldInterpolation &&
         left.horizontalAbsoluteTolerance ==
             right.horizontalAbsoluteTolerance &&
         left.verticalAbsoluteTolerance == right.verticalAbsoluteTolerance &&
         left.outputScale == right.outputScale &&
         left.outputOffset == right.outputOffset;
}

bool sameOutputFile(const WVOutputFileRecord &left,
                    const WVOutputFileRecord &right) {
  if (left.identifier != right.identifier ||
      left.destination != right.destination ||
      left.groups.size() != right.groups.size())
    return false;
  for (std::size_t index = 0; index < left.groups.size(); ++index) {
    const auto &a = left.groups[index];
    const auto &b = right.groups[index];
    if (a.identifier != b.identifier || a.name != b.name ||
        !sameSchedule(a.schedule, b.schedule) ||
        a.observerIdentifiers != b.observerIdentifiers ||
        a.containsCompleteCoefficientRestart !=
            b.containsCompleteCoefficientRestart)
      return false;
  }
  return true;
}

bool sameGraph(const WVPortableObserverRecord &left,
               const WVPortableObserverRecord &right) {
  if (left.schemaIdentifier != right.schemaIdentifier ||
      left.schemaVersion != right.schemaVersion ||
      left.stateBlocks.size() != right.stateBlocks.size() ||
      left.observers.size() != right.observers.size() ||
      left.outputFiles.size() != right.outputFiles.size())
    return false;
  for (std::size_t index = 0; index < left.stateBlocks.size(); ++index)
    if (!sameStateBlock(left.stateBlocks[index], right.stateBlocks[index]))
      return false;
  for (std::size_t index = 0; index < left.observers.size(); ++index)
    if (!sameObserver(left.observers[index], right.observers[index]))
      return false;
  for (std::size_t index = 0; index < left.outputFiles.size(); ++index)
    if (!sameOutputFile(left.outputFiles[index], right.outputFiles[index]))
      return false;
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

bool sameCursor(const WVOutputScheduleCursor &left,
                const WVOutputScheduleCursor &right) noexcept {
  return left.committedOrdinal == right.committedOrdinal &&
         sameTypedRecord(left.values, right.values);
}

bool sameObservationSchema(const WVObservationSchema &left,
                           const WVObservationSchema &right) noexcept {
  std::vector<std::uint8_t> leftBytes;
  std::vector<std::uint8_t> rightBytes;
  return encodeObservationSchemaManifest(left, leftBytes) &&
         encodeObservationSchemaManifest(right, rightBytes) &&
         leftBytes == rightBytes;
}

bool sameObservationSchemas(
    const std::vector<WVInspectedObservationSchema> &left,
    const std::vector<WVInspectedObservationSchema> &right) noexcept {
  if (left.size() != right.size())
    return false;
  for (std::size_t index = 0; index < left.size(); ++index)
    if (left[index].observerIdentifier != right[index].observerIdentifier ||
        !sameObservationSchema(left[index].schema, right[index].schema))
      return false;
  return true;
}

bool validIdentifier(const std::string &value) noexcept {
  return !value.empty() &&
         std::all_of(value.begin(), value.end(), [](unsigned char character) {
           return (character >= 'a' && character <= 'z') ||
                  (character >= 'A' && character <= 'Z') ||
                  (character >= '0' && character <= '9') || character == '-' ||
                  character == '_' || character == '.';
         });
}

bool sameTime(double first, double second) noexcept {
  const double tolerance =
      8.0 * std::numeric_limits<double>::epsilon() *
      std::max({1.0, std::abs(first), std::abs(second)});
  return std::abs(first - second) <= tolerance;
}

} // namespace

WVKernelStatus WVModelOutputGroup::evenlySpaced(
    std::string name, double outputInterval, double initialTime,
    double finalTime, WVModelOutputGroup &group, std::string identifier) {
  if (name.empty() || !std::isfinite(outputInterval) || outputInterval <= 0.0 ||
      !std::isfinite(initialTime) || std::isnan(finalTime) ||
      finalTime < initialTime)
    return invalid("An evenly spaced output group requires a name, positive "
                   "finite interval, a finite initial time, and a "
                   "nondecreasing finite or unbounded final time.");
  WVModelOutputGroup candidate;
  candidate.record_.name = std::move(name);
  candidate.record_.identifier = std::move(identifier);
  candidate.record_.schedule = {outputInterval, initialTime, finalTime};
  group = std::move(candidate);
  return WVKernelStatus::ok();
}

WVKernelStatus WVModelOutputGroup::fromRecord(WVOutputGroupRecord record,
                                              WVModelOutputGroup &group) {
  if (record.name.empty())
    return invalid("An output group record requires a name.");
  WVModelOutputGroup candidate;
  candidate.record_ = std::move(record);
  group = std::move(candidate);
  return WVKernelStatus::ok();
}

WVKernelStatus
WVModelOutputGroup::addObservingSystem(std::string observerIdentifier) {
  if (sealed_)
    return invalid("The output group is sealed.");
  if (observerIdentifier.empty())
    return invalid("An observing-system identifier is required.");
  if (std::find(record_.observerIdentifiers.begin(),
                record_.observerIdentifiers.end(), observerIdentifier) !=
      record_.observerIdentifiers.end())
    return invalid("The observing system is already a member of this group.");
  record_.observerIdentifiers.push_back(std::move(observerIdentifier));
  return WVKernelStatus::ok();
}

WVKernelStatus
WVModelOutputGroup::containsCompleteCoefficientRestart(bool value) {
  if (sealed_)
    return invalid("The output group is sealed.");
  record_.containsCompleteCoefficientRestart = value;
  return WVKernelStatus::ok();
}

WVKernelStatus WVModelOutputGroup::scheduleContinuation(
    WVOutputScheduleCursor continuation) {
  if (sealed_)
    return invalid("The output group is sealed.");
  if (continuation.committedOrdinal < WVNoCommittedOutputOrdinal)
    return invalid("An output schedule ordinal cannot be less than -1.");
  continuation_ = std::move(continuation);
  hasExplicitContinuation_ = true;
  return WVKernelStatus::ok();
}

std::size_t WVModelOutputGroup::persistentBytes() const noexcept {
  std::size_t result = sizeof(*this) + record_.identifier.capacity() +
                       record_.name.capacity() +
                       record_.observerIdentifiers.capacity() *
                           sizeof(std::string);
  for (const auto &identifier : record_.observerIdentifiers)
    result += identifier.capacity();
  return result;
}

WVKernelStatus WVModelOutputFile::create(std::string destination,
                                         WVModelOutputFile &file,
                                         std::string identifier) {
  if (destination.empty())
    return invalid("An output destination is required.");
  WVModelOutputFile candidate;
  try {
    candidate.destination_ = normalizedDestination(destination);
  } catch (const std::exception &error) {
    return invalid(std::string("The output destination cannot be resolved: ") +
                   error.what());
  }
  candidate.identifier_ =
      identifier.empty()
          ? stableIdentifier("output-file", candidate.destination_)
          : std::move(identifier);
  file = std::move(candidate);
  return WVKernelStatus::ok();
}

WVKernelStatus WVModelOutputFile::addOutputGroup(WVModelOutputGroup group) {
  if (sealed_)
    return invalid("The output file is sealed.");
  if (group.sealed_)
    return invalid("A sealed output group cannot be added.");
  if (group.record_.identifier.empty())
    group.record_.identifier =
        stableIdentifier("output-group", identifier_ + ":" + group.name());
  if (outputGroupWithName(group.name()) != nullptr)
    return invalid("The output file already contains a group with this name.");
  const auto duplicateIdentifier =
      std::find_if(groups_.begin(), groups_.end(), [&](const auto &existing) {
        return existing.identifier() == group.identifier();
      });
  if (duplicateIdentifier != groups_.end())
    return invalid("The output file already contains this group identifier.");
  groups_.push_back(std::move(group));
  return WVKernelStatus::ok();
}

WVKernelStatus WVModelOutputFile::addNewEvenlySpacedOutputGroup(
    std::string name, double outputInterval, double initialTime,
    double finalTime, WVModelOutputGroup *&group, std::string identifier) {
  group = nullptr;
  WVModelOutputGroup candidate;
  auto status = WVModelOutputGroup::evenlySpaced(
      std::move(name), outputInterval, initialTime, finalTime, candidate,
      std::move(identifier));
  if (!status)
    return status;
  status = addOutputGroup(std::move(candidate));
  if (!status)
    return status;
  group = &groups_.back();
  return WVKernelStatus::ok();
}

WVKernelStatus
WVModelOutputFile::addObservingSystem(std::string observerIdentifier) {
  if (sealed_)
    return invalid("The output file is sealed.");
  if (groups_.size() != 1)
    return invalid("File-level observing-system addition requires exactly "
                   "one output group.");
  return groups_.front().addObservingSystem(std::move(observerIdentifier));
}

WVModelOutputGroup *
WVModelOutputFile::outputGroupWithName(const std::string &name) noexcept {
  const auto found = std::find_if(groups_.begin(), groups_.end(),
                                  [&](const auto &group) {
                                    return group.name() == name;
                                  });
  return found == groups_.end() ? nullptr : &*found;
}

const WVModelOutputGroup *WVModelOutputFile::outputGroupWithName(
    const std::string &name) const noexcept {
  const auto found = std::find_if(groups_.begin(), groups_.end(),
                                  [&](const auto &group) {
                                    return group.name() == name;
                                  });
  return found == groups_.end() ? nullptr : &*found;
}

WVModelOutputGroup *WVModelOutputFile::outputGroupWithIdentifier(
    const std::string &identifier) noexcept {
  const auto found = std::find_if(groups_.begin(), groups_.end(),
                                  [&](const auto &group) {
                                    return group.identifier() == identifier;
                                  });
  return found == groups_.end() ? nullptr : &*found;
}

const WVModelOutputGroup *WVModelOutputFile::outputGroupWithIdentifier(
    const std::string &identifier) const noexcept {
  const auto found = std::find_if(groups_.begin(), groups_.end(),
                                  [&](const auto &group) {
                                    return group.identifier() == identifier;
                                  });
  return found == groups_.end() ? nullptr : &*found;
}

std::size_t WVModelOutputFile::persistentBytes() const noexcept {
  std::size_t result = sizeof(*this) + destination_.capacity() +
                       identifier_.capacity();
  for (const auto &group : groups_)
    result += group.persistentBytes();
  return result;
}

class WVModelOutputConfiguration::Impl final {
public:
  std::shared_ptr<const WVExtensionCatalog> catalog;
  WVPortableObserverDescriptor descriptor;
  WVOutputPlan plan;
  std::vector<WVOutputScheduleContinuation> scheduleContinuations;
  std::vector<WVOutputDestinationProgress> destinationProgress;
  std::vector<WVInspectedObservationSchema> observationSchemas;
  WVModelOutputPolicy policy = WVModelOutputPolicy::create;
};

WVModelOutputConfiguration::WVModelOutputConfiguration()
    : impl_(std::make_unique<Impl>()) {}
WVModelOutputConfiguration::~WVModelOutputConfiguration() = default;
WVModelOutputConfiguration::WVModelOutputConfiguration(
    WVModelOutputConfiguration &&) noexcept = default;
WVModelOutputConfiguration &WVModelOutputConfiguration::operator=(
    WVModelOutputConfiguration &&) noexcept = default;

WVKernelStatus WVModelOutputConfiguration::build(
    WVPortableObserverRecord observerRecord,
    std::vector<WVModelOutputFile> files, WVModelOutputPolicy policy,
    std::shared_ptr<const WVExtensionCatalog> catalog,
    double initialTime, double finalTime,
    WVModelOutputConfiguration &configuration) {
  if (!catalog)
    return invalid("Output configuration requires an extension catalog.");
  if (files.empty())
    return invalid("At least one output file is required.");
  try {
    observerRecord.outputFiles.clear();
    std::vector<WVOutputScheduleContinuation> continuations;
    bool hasAnyExplicitContinuation = false;
    bool hasAnyImplicitContinuation = false;
    for (auto &file : files) {
      if (file.sealed_)
        return invalid("An output file builder was already consumed.");
      file.destination_ = normalizedDestination(file.destination_);
      WVOutputFileRecord record{file.identifier_, file.destination_, {}};
      for (auto &group : file.groups_) {
        if (group.sealed_)
          return invalid("An output group builder was already consumed.");
        if (group.record_.identifier.empty())
          group.record_.identifier = stableIdentifier(
              "output-group", file.identifier_ + ":" + group.name());
        record.groups.push_back(group.record_);
        continuations.push_back({file.identifier_, group.identifier(),
                                 group.continuation_});
        hasAnyExplicitContinuation |= group.hasExplicitContinuation_;
        hasAnyImplicitContinuation |= !group.hasExplicitContinuation_;
        group.sealed_ = true;
      }
      file.sealed_ = true;
      observerRecord.outputFiles.push_back(std::move(record));
    }
    if (hasAnyExplicitContinuation && hasAnyImplicitContinuation)
      return invalid("Schedule continuation must be supplied for every output "
                     "group or for none of them.");
    files.clear();
    if (!hasAnyExplicitContinuation)
      continuations.clear();
    return compile(std::move(observerRecord), {}, std::move(continuations),
                   policy, std::move(catalog), initialTime, finalTime,
                   configuration);
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Output configuration allocation failed."};
  } catch (const std::exception &error) {
    return invalid(std::string("Output configuration failed: ") + error.what());
  }
}

WVKernelStatus WVModelOutputConfiguration::compile(
    WVPortableObserverRecord observerRecord,
    std::vector<WVInspectedObservationSchema> observationSchemas,
    std::vector<WVOutputScheduleContinuation> scheduleContinuations,
    WVModelOutputPolicy policy,
    std::shared_ptr<const WVExtensionCatalog> catalog, double initialTime,
    double finalTime, WVModelOutputConfiguration &configuration,
    const WVTransformConstantStratificationConfiguration
        *planningConfiguration,
    bool isDynamicsLinear,
    const WVTransformStateDescription *planningStateDescription) {
  if (!catalog)
    return invalid("Output configuration requires an extension catalog.");
  if (!std::isfinite(initialTime) || !std::isfinite(finalTime) ||
      finalTime < initialTime)
    return invalid("Output compilation requires a finite nondecreasing time "
                   "interval.");
  if (observerRecord.outputFiles.empty())
    return invalid("At least one output file is required.");
  try {
    if (observerRecord.schemaIdentifier !=
            WVPortableObserverContractIdentifier ||
        observerRecord.schemaVersion != WVPortableObserverContractVersion)
      return invalid("Unsupported portable observing-system contract schema.");

    std::set<std::string> blockIdentifiers;
    const WVStateBlockRecord *canonicalAp = nullptr;
    for (const auto &block : observerRecord.stateBlocks) {
      if (!validIdentifier(block.identifier) ||
          !blockIdentifiers.insert(block.identifier).second)
        return invalid("State-block identifiers must be unique portable "
                       "identifiers.");
      if (block.identifier == "Ap")
        canonicalAp = &block;
    }
    WVIntegrationStateLayout rawLayout;
    WVKernelStatus status;
    if (planningStateDescription == nullptr) {
      if (canonicalAp == nullptr || canonicalAp->dimensions.size() != 2)
        return invalid("Output compilation requires a canonical [Nj,Nkl] Ap "
                       "state block.");
      status = WVIntegrationStateLayout::create(
          {canonicalAp->dimensions[0], canonicalAp->dimensions[1]},
          observerRecord, rawLayout);
    } else {
      status = WVIntegrationStateLayout::create(
          *planningStateDescription, observerRecord, rawLayout);
    }
    if (!status)
      return status;

    std::set<std::string> observerIdentifiers;
    for (auto &observer : observerRecord.observers) {
      if (!validIdentifier(observer.identifier) || observer.name.empty() ||
          !observerIdentifiers.insert(observer.identifier).second)
        return invalid("Observer identifiers must be unique portable "
                       "identifiers and observer names must be nonempty.");
      if (catalog->observers().registration(observer.typeIdentifier,
                                            observer.contractVersion) ==
          nullptr)
        return {WVKernelStatusCode::unsupportedOperation,
                "Unsupported observing-system identity or contract version."};
      WVPortableTypedRecord canonicalConfiguration;
      status = catalog->observers().resolveConfiguration(
          observer, canonicalConfiguration);
      if (!status)
        return status;
      observer.configuration = std::move(canonicalConfiguration);
      std::set<std::string> referencedBlocks;
      for (const auto &identifier : observer.stateBlockIdentifiers)
        if (blockIdentifiers.find(identifier) == blockIdentifiers.end() ||
            !referencedBlocks.insert(identifier).second)
          return invalid("Observer state-block references must be unique and "
                         "resolve within the canonical graph.");
    }

    std::set<std::string> declaredSchemaObservers;
    for (const auto &declared : observationSchemas) {
      if (observerIdentifiers.find(declared.observerIdentifier) ==
              observerIdentifiers.end() ||
          !declaredSchemaObservers.insert(declared.observerIdentifier).second)
        return invalid("Declared observation schemas must reference distinct "
                       "canonical observers.");
      status = validateObservationSchema(declared.schema);
      if (!status)
        return status;
    }
    if (!observationSchemas.empty()) {
      WVObserverOutputPlanningContext planningContext;
      planningContext.configuration = planningConfiguration;
      planningContext.stateBlocks = observerRecord.stateBlocks.data();
      planningContext.stateBlockCount = observerRecord.stateBlocks.size();
      planningContext.isDynamicsLinear = isDynamicsLinear;
      planningContext.stateLayout = &rawLayout;
      for (const auto &observer : observerRecord.observers) {
        const auto declared = std::find_if(
            observationSchemas.begin(), observationSchemas.end(),
            [&](const auto &candidate) {
              return candidate.observerIdentifier == observer.identifier;
            });
        if (declared == observationSchemas.end())
          continue;
        WVObserverOutputPlan plan;
        status = catalog->observers().resolveOutputPlan(
            observer, planningContext, plan);
        if (!status)
          return status;
        status = validateObservationSchema(plan.schema);
        if (!status)
          return status;
        if (!sameObservationSchema(declared->schema, plan.schema))
          return invalid("The registered observer output plan differs from "
                         "the canonical schema restored from NetCDF.");
      }
    }

    std::set<std::string> destinations;
    std::set<std::string> fileIdentifiers;
    std::size_t groupCount = 0;
    std::vector<std::string> destinationPaths;
    for (auto &file : observerRecord.outputFiles) {
      file.destination = normalizedDestination(file.destination);
      if (!validIdentifier(file.identifier) ||
          !fileIdentifiers.insert(file.identifier).second ||
          !destinations.insert(file.destination).second)
        return invalid("Output file identifiers and destinations must be "
                       "unique and nonempty.");
      std::error_code error;
      const bool exists = std::filesystem::exists(file.destination, error);
      if (error)
        return invalid("An output destination cannot be inspected: " +
                       error.message());
      if (policy == WVModelOutputPolicy::create && exists)
        return invalid("Create policy requires every destination to be absent.");
      if (policy != WVModelOutputPolicy::create && !exists)
        return invalid("Replace and append policies require every destination "
                       "to exist.");
      destinationPaths.push_back(file.destination);
      std::set<std::string> groupIdentifiers;
      std::set<std::string> groupNames;
      std::size_t restartGroupCount = 0;
      for (const auto &group : file.groups) {
        if (!validIdentifier(group.identifier) || group.name.empty() ||
            !groupIdentifiers.insert(group.identifier).second ||
            !groupNames.insert(group.name).second)
          return invalid("Output group identifiers and names must be unique "
                         "within each file.");
        std::set<std::string> groupObservers;
        for (const auto &identifier : group.observerIdentifiers)
          if (observerIdentifiers.find(identifier) == observerIdentifiers.end() ||
              !groupObservers.insert(identifier).second)
            return invalid("Output group observer membership must be unique "
                           "and resolve within the canonical graph.");
        restartGroupCount +=
            group.containsCompleteCoefficientRestart ? 1U : 0U;
        ++groupCount;
      }
      if (restartGroupCount != 1)
        return invalid("Every output file must designate exactly one complete "
                       "coefficient-restart group.");
    }

    WVModelOutputNetCDFInspection appendInspection;
    std::vector<WVOutputDestinationProgress> destinationProgress;
    if (policy == WVModelOutputPolicy::append) {
      const auto inspected = WVModelOutputNetCDFSink::inspect(
          destinationPaths, *catalog, appendInspection);
      if (!inspected)
        return fromCheckpoint(inspected);
      for (auto &observer : appendInspection.observerRecord.observers) {
        WVPortableTypedRecord canonicalConfiguration;
        status = catalog->observers().resolveConfiguration(
            observer, canonicalConfiguration);
        if (!status)
          return status;
        observer.configuration = std::move(canonicalConfiguration);
      }
      if (!sameGraph(observerRecord, appendInspection.observerRecord))
        return invalid("Append destinations do not contain the complete "
                       "configured observer and output graph.");
      if (observationSchemas.empty())
        observationSchemas = appendInspection.observationSchemas;
      else if (!sameObservationSchemas(observationSchemas,
                                       appendInspection.observationSchemas))
        return invalid("Append destinations contain incompatible declared "
                       "observation schemas.");
      destinationProgress = appendInspection.destinationProgress;
    } else {
      destinationProgress.reserve(groupCount);
      for (const auto &file : observerRecord.outputFiles)
        for (const auto &group : file.groups) {
          WVOutputDestinationProgress progress;
          progress.fileIdentifier = file.identifier;
          progress.groupIdentifier = group.identifier;
          destinationProgress.push_back(std::move(progress));
        }
    }

    if (scheduleContinuations.empty()) {
      scheduleContinuations.reserve(groupCount);
      std::size_t destinationIndex = 0;
      for (const auto &file : observerRecord.outputFiles)
        for (const auto &group : file.groups) {
          WVOutputScheduleCursor cursor;
          if (policy == WVModelOutputPolicy::append)
            cursor = destinationProgress[destinationIndex]
                         .committedScheduleCursor;
          scheduleContinuations.push_back(
              {file.identifier, group.identifier, std::move(cursor)});
          ++destinationIndex;
        }
    }
    if (scheduleContinuations.size() != groupCount ||
        destinationProgress.size() != groupCount)
      return invalid("Schedule continuation and destination progress must each "
                     "contain one entry per output group.");

    std::size_t groupIndex = 0;
    for (const auto &file : observerRecord.outputFiles) {
      for (const auto &group : file.groups) {
        const auto &continuation = scheduleContinuations[groupIndex];
        const auto &progress = destinationProgress[groupIndex];
        if (continuation.fileIdentifier != file.identifier ||
            continuation.groupIdentifier != group.identifier ||
            progress.fileIdentifier != file.identifier ||
            progress.groupIdentifier != group.identifier)
          return invalid("Schedule continuation or destination progress does "
                         "not match deterministic file/group order.");
        std::shared_ptr<const WVOutputSchedule> schedule;
        status = catalog->outputSchedules().resolve(group.schedule, schedule);
        if (!status)
          return status;
        status = schedule->validateCursor(continuation.cursor);
        if (!status)
          return status;
        double continuationCommittedTime = 0.0;
        bool hasContinuationCommittedTime = false;
        status = schedule->committedTime(continuation.cursor,
                                         continuationCommittedTime,
                                         hasContinuationCommittedTime);
        if (!status)
          return status;
        if (hasContinuationCommittedTime &&
            continuationCommittedTime > initialTime &&
            !sameTime(continuationCommittedTime, initialTime))
          return invalid("Schedule continuation lies beyond the selected "
                         "restart time.");
        if (policy == WVModelOutputPolicy::append) {
          status = schedule->validateCursor(progress.committedScheduleCursor);
          if (!status)
            return status;
          double destinationCommittedTime = 0.0;
          bool hasDestinationCommittedTime = false;
          status = schedule->committedTime(progress.committedScheduleCursor,
                                           destinationCommittedTime,
                                           hasDestinationCommittedTime);
          if (!status)
            return status;
          if (!sameCursor(continuation.cursor,
                          progress.committedScheduleCursor))
            return invalid("Append schedule continuation cursor differs from "
                           "the destination's complete typed cursor (source " +
                           std::to_string(
                               continuation.cursor.committedOrdinal) +
                           ", destination " +
                           std::to_string(progress.committedScheduleCursor
                                              .committedOrdinal) +
                           ").");
          if (hasDestinationCommittedTime != progress.hasCommittedTime)
            return invalid("Append schedule state and destination time-last "
                           "marker availability differ.");
          if (progress.hasCommittedTime &&
              !sameTime(destinationCommittedTime,
                        progress.lastCommittedTime))
            return invalid("Append schedule time differs from the destination "
                           "time-last commit marker.");
          if ((progress.hasCommittedTime && progress.recordCount == 0) ||
              (!progress.hasCommittedTime && progress.recordCount != 0))
            return invalid("Append destination record count and time-last "
                           "commit marker disagree.");
          for (const auto &axis : progress.unlimitedAxes)
            if (axis.committedCount != axis.physicalCount)
              return invalid("Append ragged-axis committed and physical "
                             "offsets differ.");
        }
        ++groupIndex;
      }
    }

    // Runtime implementations are constructed only after every raw graph,
    // cursor, destination, and append-progress check above has passed.
    auto candidate = std::make_unique<Impl>();
    candidate->catalog = std::move(catalog);
    status = WVPortableObserverDescriptor::create(
        observerRecord, candidate->catalog, candidate->descriptor);
    if (!status)
      return status;
    WVIntegrationStateLayout resolvedLayout;
    if (planningStateDescription == nullptr) {
      status = WVIntegrationStateLayout::create(
          rawLayout.coefficientShape(), candidate->descriptor,
          resolvedLayout);
    } else {
      status = WVIntegrationStateLayout::create(
          *planningStateDescription, candidate->descriptor, resolvedLayout);
    }
    if (!status)
      return status;
    status = WVOutputPlan::create(
        resolvedLayout, candidate->descriptor, candidate->catalog,
        initialTime, finalTime, scheduleContinuations, candidate->plan);
    if (!status)
      return status;
    candidate->scheduleContinuations =
        candidate->plan.initialContinuations();
    candidate->destinationProgress = std::move(destinationProgress);
    candidate->observationSchemas = std::move(observationSchemas);
    candidate->policy = policy;
    configuration.impl_ = std::move(candidate);
    return WVKernelStatus::ok();
  } catch (const std::bad_alloc &) {
    return {WVKernelStatusCode::allocationFailure,
            "Output configuration allocation failed."};
  } catch (const std::exception &error) {
    return invalid(std::string("Output configuration failed: ") + error.what());
  }
}

WVCheckpointStatus WVModelOutputConfiguration::openNetCDFSink(
    const WVModelOutputNetCDFConfiguration &configuration,
    const WVIntegrationStateLayout &stateLayout,
    WVObserverSampleSource *sampleSource, WVModelOutputNetCDFSink &sink) const {
  switch (impl_->policy) {
  case WVModelOutputPolicy::create:
    return WVModelOutputNetCDFSink::createNew(configuration, impl_->descriptor,
                                              impl_->plan, stateLayout,
                                              sampleSource, sink);
  case WVModelOutputPolicy::replace:
    return WVModelOutputNetCDFSink::replaceExisting(
        configuration, impl_->descriptor, impl_->plan, stateLayout,
        sampleSource, sink);
  case WVModelOutputPolicy::append:
    return WVModelOutputNetCDFSink::openAppend(
        configuration, impl_->descriptor, impl_->plan, stateLayout,
        sampleSource, impl_->destinationProgress, sink);
  }
  return {
      WVCheckpointStatusCode::invalidValue, "Unsupported output policy.", {}};
}

const WVPortableObserverDescriptor &
WVModelOutputConfiguration::descriptor() const noexcept {
  return impl_->descriptor;
}

const std::shared_ptr<const WVExtensionCatalog> &
WVModelOutputConfiguration::catalog() const noexcept {
  return impl_->catalog;
}

const WVOutputPlan &WVModelOutputConfiguration::plan() const noexcept {
  return impl_->plan;
}

const std::vector<WVOutputScheduleContinuation> &
WVModelOutputConfiguration::scheduleContinuations() const noexcept {
  return impl_->scheduleContinuations;
}

const std::vector<WVOutputDestinationProgress> &
WVModelOutputConfiguration::destinationProgress() const noexcept {
  return impl_->destinationProgress;
}

const std::vector<WVInspectedObservationSchema> &
WVModelOutputConfiguration::observationSchemas() const noexcept {
  return impl_->observationSchemas;
}

WVModelOutputPolicy WVModelOutputConfiguration::policy() const noexcept {
  return impl_->policy;
}

std::size_t WVModelOutputConfiguration::persistentBytes() const noexcept {
  // The plan shares and accounts the same immutable descriptor graph; the
  // configuration's descriptor wrapper itself is already part of sizeof(Impl).
  std::size_t bytes = sizeof(*this) + sizeof(Impl) +
                      impl_->plan.persistentBytes() - sizeof(impl_->plan) +
                      impl_->scheduleContinuations.capacity() *
                          sizeof(WVOutputScheduleContinuation) +
                      impl_->destinationProgress.capacity() *
                          sizeof(WVOutputDestinationProgress) +
                      impl_->observationSchemas.capacity() *
                          sizeof(WVInspectedObservationSchema);
  for (const auto &continuation : impl_->scheduleContinuations)
    bytes += continuation.fileIdentifier.capacity() +
             continuation.groupIdentifier.capacity() +
             continuation.cursor.values.persistentBytes() -
                 sizeof(WVPortableTypedRecord);
  for (const auto &progress : impl_->destinationProgress) {
    bytes += progress.fileIdentifier.capacity() +
             progress.groupIdentifier.capacity() +
             progress.committedScheduleCursor.values.persistentBytes() -
                 sizeof(WVPortableTypedRecord) +
             progress.unlimitedAxes.capacity() *
                 sizeof(WVOutputDestinationAxisProgress);
    for (const auto &axis : progress.unlimitedAxes)
      bytes += axis.axisIdentifier.capacity();
  }
  for (const auto &schema : impl_->observationSchemas)
    bytes += schema.observerIdentifier.capacity() +
             observationSchemaRetainedBytes(schema.schema);
  return bytes;
}

} // namespace wavevortex::runtime
