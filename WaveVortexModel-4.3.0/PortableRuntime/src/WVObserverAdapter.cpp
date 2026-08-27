#include "WVObserverAdapter.hpp"
#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVObserverOutputProvider.hpp"
#include "WVLegacyObserverCompatibility.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <tuple>
#include <utility>

namespace wavevortex::runtime::detail {
namespace {

WVKernelStatus invalid(std::string message) {
  return {WVKernelStatusCode::invalidConfiguration, std::move(message)};
}

WVKernelStatus configuredText(const WVObserverRecord &observer,
                              const char *name,
                              const std::vector<std::string> &fallback,
                              std::vector<std::string> &output) {
  const auto *value = observer.configuration.value(name);
  if (value == nullptr) {
    output = fallback;
    return WVKernelStatus::ok();
  }
  const auto *text =
      std::get_if<std::vector<std::string>>(&value->storage);
  if (text == nullptr)
    return invalid(std::string("Observer configuration value '") + name +
                   "' must contain text.");
  output = *text;
  return WVKernelStatus::ok();
}

using LegacyConfigurationPopulator =
    WVKernelStatus (*)(const WVObserverRecord &, WVPortableTypedRecord &);

void addLegacyText(WVPortableTypedRecord &record, std::string name,
                   const std::vector<std::string> &values) {
  if (!values.empty())
    record.values.push_back({std::move(name), {values.size()}, values});
}

void addLegacyReal(WVPortableTypedRecord &record, std::string name,
                   const std::vector<double> &values) {
  if (!values.empty())
    record.values.push_back({std::move(name), {values.size()}, values});
}

void addLegacyBoolean(WVPortableTypedRecord &record, std::string name,
                      bool value) {
  record.values.push_back(
      {std::move(name), {},
       std::vector<std::uint8_t>{static_cast<std::uint8_t>(value ? 1 : 0)}});
}

void addLegacyScalar(WVPortableTypedRecord &record, std::string name,
                     double value) {
  record.values.push_back(
      {std::move(name), {}, std::vector<double>{value}});
}

WVKernelStatus resolveLegacyConfiguration(
    const WVObserverRecord &observer, LegacyConfigurationPopulator populate,
    WVPortableTypedRecord &configuration) {
  if (!observer.configuration.schemaIdentifier.empty()) {
    const auto status = validatePortableTypedRecord(
        observer.configuration, {1024 * 1024, true, true});
    if (!status)
      return status;
    configuration = observer.configuration;
    return WVKernelStatus::ok();
  }
  WVPortableTypedRecord candidate;
  candidate.schemaIdentifier =
      "legacy-" + observer.typeIdentifier + "-configuration-v1";
  candidate.schemaVersion = 1;
  const auto populateStatus = populate(observer, candidate);
  if (!populateStatus)
    return populateStatus;
  const auto status =
      validatePortableTypedRecord(candidate, {1024 * 1024, true, true});
  if (!status)
    return status;
  configuration = std::move(candidate);
  return WVKernelStatus::ok();
}

WVKernelStatus populateCoefficientConfiguration(
    const WVObserverRecord &, WVPortableTypedRecord &) {
  return WVKernelStatus::ok();
}

WVKernelStatus populateEulerianConfiguration(
    const WVObserverRecord &observer, WVPortableTypedRecord &record) {
  addLegacyText(record, "fieldNames", observer.fieldNames);
  return WVKernelStatus::ok();
}

WVKernelStatus populateMooringConfiguration(
    const WVObserverRecord &observer, WVPortableTypedRecord &record) {
  addLegacyText(record, "trackedFieldNames", observer.fieldNames);
  addLegacyReal(record, "x", observer.x);
  addLegacyReal(record, "y", observer.y);
  addLegacyReal(record, "z", observer.z);
  return WVKernelStatus::ok();
}

WVKernelStatus populateParticleConfiguration(
    const WVObserverRecord &observer, WVPortableTypedRecord &record) {
  addLegacyText(record, "trackedFieldNames", observer.fieldNames);
  addLegacyReal(record, "x", observer.x);
  addLegacyReal(record, "y", observer.y);
  addLegacyReal(record, "z", observer.z);
  addLegacyBoolean(record, "isXYOnly", observer.isXYOnly);
  const auto interpolation = [](WVPositionInterpolation value) {
    return value == WVPositionInterpolation::spline ? std::string("spline")
                                                     : std::string("linear");
  };
  record.values.push_back(
      {"advectionInterpolation", {},
       std::vector<std::string>{
           interpolation(observer.advectionInterpolation)}});
  record.values.push_back(
      {"trackedFieldInterpolation", {},
       std::vector<std::string>{
           interpolation(observer.trackedFieldInterpolation)}});
  addLegacyScalar(record, "horizontalAbsoluteTolerance",
                  observer.horizontalAbsoluteTolerance);
  addLegacyScalar(record, "verticalAbsoluteTolerance",
                  observer.verticalAbsoluteTolerance);
  return WVKernelStatus::ok();
}

WVKernelStatus populateTracerConfiguration(
    const WVObserverRecord &observer, WVPortableTypedRecord &record) {
  addLegacyBoolean(record, "isXYOnly", observer.isXYOnly);
  addLegacyBoolean(record, "shouldAntialias", observer.shouldAntialias);
  return WVKernelStatus::ok();
}

WVKernelStatus invokeLegacyOperation(
    const std::function<WVKernelStatus()> &operation,
    const char *description) {
  return operation
             ? operation()
             : invalid(std::string("The runtime cannot bind the legacy ") +
                       description + " observer operation.");
}

WVKernelStatus selectFullFieldOperation(
    const WVObserverRecord &,
    const WVLegacyObserverOperationBinder &binder) {
  return invokeLegacyOperation(binder.fullField, "full-field");
}

WVKernelStatus selectFixedVerticalProfileOperation(
    const WVObserverRecord &,
    const WVLegacyObserverOperationBinder &binder) {
  return invokeLegacyOperation(binder.fixedVerticalProfiles,
                               "fixed-profile");
}

WVKernelStatus selectMovingPositionOperation(
    const WVObserverRecord &,
    const WVLegacyObserverOperationBinder &binder) {
  return invokeLegacyOperation(binder.movingPositions, "moving-position");
}

WVKernelStatus selectIntegratedStateOperation(
    const WVObserverRecord &,
    const WVLegacyObserverOperationBinder &binder) {
  return invokeLegacyOperation(binder.integratedState, "integrated-state");
}

void addAxis(WVObservationSchema &schema, const std::string &name,
             std::size_t extent, WVObservationCoordinateRole role) {
  const auto found =
      std::find_if(schema.axes.begin(), schema.axes.end(),
                   [&](const auto &axis) { return axis.identifier == name; });
  if (found == schema.axes.end())
    schema.axes.push_back(
        {name, name, WVObservationAxisKind::fixed, extent, role});
}

WVObservationMetadataVariable metadataReal(const std::string &name,
                                           double value) {
  WVObservationMetadataVariable variable;
  variable.name = name;
  variable.value = WVObservationValue::ownReal(name, {}, {value});
  return variable;
}

WVObservationMetadataVariable metadataBoolean(const std::string &name,
                                              bool value) {
  WVObservationMetadataVariable variable;
  variable.name = name;
  variable.value = WVObservationValue::ownBoolean(
      name, {}, {static_cast<std::uint8_t>(value ? 1 : 0)});
  variable.isLogicalType = true;
  return variable;
}

void addConstantReal(
    WVObserverOutputPlan &plan, std::string identifier, std::string name,
    std::vector<std::string> dimensions, WVObservationValueLayout layout,
    std::vector<double> values, std::string units, std::string description,
    WVObservationCoordinateRole role,
    std::vector<WVObservationAttribute> attributes = {}) {
  WVObservationVariable variable;
  variable.identifier = identifier;
  variable.name = std::move(name);
  variable.dimensionIdentifiers = std::move(dimensions);
  variable.layout = layout;
  variable.units = std::move(units);
  variable.description = std::move(description);
  variable.attributes = std::move(attributes);
  variable.coordinateRole = role;
  plan.schema.variables.push_back(std::move(variable));
  std::vector<std::size_t> extents;
  for (const auto &dimension :
       plan.schema.variables.back().dimensionIdentifiers) {
    const auto axis = std::find_if(
        plan.schema.axes.begin(), plan.schema.axes.end(),
        [&](const auto &candidate) {
          return candidate.identifier == dimension;
        });
    extents.push_back(axis == plan.schema.axes.end() ? 0 : axis->extent);
  }
  plan.constantValues.push_back(WVObservationValue::ownReal(
      std::move(identifier), std::move(extents), std::move(values)));
}

void addChannel(WVObserverOutputPlan &plan, WVObservationVariable variable,
                WVObserverOutputChannel channel) {
  channel.variableIdentifier = variable.identifier;
  plan.schema.variables.push_back(std::move(variable));
  plan.channels.push_back(std::move(channel));
}

WVObservationVariable fieldVariable(
    const WVPortableVariableMetadata &metadata, std::string identifier,
    std::string name, std::vector<std::string> dimensions,
    WVObservationValueLayout layout, std::string descriptionSuffix = {}) {
  WVObservationVariable variable;
  variable.identifier = std::move(identifier);
  variable.name = std::move(name);
  variable.scalarType = metadata.kind == WVPortableVariableKind::coefficient
                            ? WVObservationScalarType::complex64
                            : WVObservationScalarType::real64;
  variable.dimensionIdentifiers = std::move(dimensions);
  variable.layout = layout;
  variable.units = metadata.units;
  variable.description = std::string(metadata.description) + descriptionSuffix;
  if (metadata.netCDFAttributeCount != 0)
    variable.attributes.push_back(
        {metadata.netCDFAttribute.name, metadata.netCDFAttribute.value});
  return variable;
}

WVKernelStatus buildLegacyOutputPlan(
    const std::string &typeIdentifier,
    const WVObserverRecord &observer,
    const WVObserverExecutionPlan &execution,
    const WVLegacyObserverOperationResolver &resolveOperation,
    const WVObserverOutputPlanningContext &context,
    WVObserverOutputPlan &output) {
  if (context.stateLayout == nullptr)
    return invalid("Observer output planning requires a resolved state layout.");
  const auto *configuration = context.configuration;
  const auto &spatialDimensions = context.stateLayout->spatialDimensions();
  if (spatialDimensions.size() < 2)
    return invalid("Observer output planning requires two horizontal dimensions.");
  WVObserverOutputPlan plan;
  plan.schema.identifier =
      "legacy-" + observer.identifier + "-observation-v1";
  plan.schema.preservesLegacyEncoding = true;
  plan.schema.metadata.attributes.push_back(
      {"AnnotatedClass", typeIdentifier});
  plan.schema.metadata.attributes.push_back(
      {"portableIdentifier", observer.identifier});
  if (!execution.persistedName.empty())
    plan.schema.metadata.attributes.push_back(
        {"name", execution.persistedName});
  if (!execution.fieldListAttribute.empty())
    plan.schema.metadata.stringListAttributes.push_back(
        {execution.fieldListAttribute, execution.outputFields});

  const auto outputLayout = [&](const WVPortableVariableMetadata &metadata) {
    return context.isDynamicsLinear &&
                   !metadata.isVariableWithLinearTimeStep
               ? WVObservationValueLayout::initialValue
               : WVObservationValueLayout::record;
  };
  const auto addFullField = [&](const std::string &field) {
    const auto *metadata = findPortableVariable(field);
    if (metadata == nullptr)
      return invalid("Unsupported observer field: " + field + ".");
    std::vector<std::string> names;
    std::vector<std::size_t> dimensions;
    WVObserverOutputChannel channel;
    channel.sourceIdentifier = field;
    if (metadata->kind == WVPortableVariableKind::coefficient) {
      channel.source = WVObserverOutputChannelSource::coefficient;
      const auto family = std::find_if(
          context.stateLayout->coefficientFamilies().begin(),
          context.stateLayout->coefficientFamilies().end(),
          [&](const auto &candidate) { return candidate.identifier == field; });
      if (family == context.stateLayout->coefficientFamilies().end())
        return invalid("Coefficient field is absent from the resolved transform: " +
                       field + ".");
      channel.coefficientFamily = static_cast<std::size_t>(
          family - context.stateLayout->coefficientFamilies().begin());
      dimensions = family->spectralDimensions;
      names.reserve(dimensions.size());
      if (dimensions.size() == 1)
        names.emplace_back("kl");
      else if (dimensions.size() == 2) {
        names.emplace_back("j");
        names.emplace_back("kl");
      } else
        for (std::size_t dimension = 0; dimension < dimensions.size();
             ++dimension)
          names.emplace_back("coefficient_dimension_" +
                             std::to_string(dimension + 1));
    } else if (metadata->kind == WVPortableVariableKind::field) {
      if (metadata->naturalRank == WVPortableNaturalRank::scalar) {
        dimensions = {};
      } else if (spatialDimensions.size() == 2) {
        dimensions = spatialDimensions;
        names = {"x", "y"};
      } else {
        for (std::size_t dimension = 0; dimension < metadata->dimensionCount;
             ++dimension)
          names.emplace_back(metadata->dimensions[dimension]);
        if (metadata->naturalRank == WVPortableNaturalRank::vertical)
          dimensions = {spatialDimensions[2]};
        else if (metadata->naturalRank == WVPortableNaturalRank::horizontal)
          dimensions = {spatialDimensions[0], spatialDimensions[1]};
        else
          dimensions = spatialDimensions;
      }
      channel.source = WVObserverOutputChannelSource::sampledField;
    } else {
      return invalid("Unsupported observer field: " + field + ".");
    }
    for (std::size_t index = 0; index < names.size(); ++index)
      addAxis(plan.schema, names[index], dimensions[index],
              WVObservationCoordinateRole::none);
    addChannel(plan,
               fieldVariable(*metadata, "derived-" + field, field, names,
                             outputLayout(*metadata)),
               std::move(channel));
    return WVKernelStatus::ok();
  };

  const bool isCompleteCoefficientObserver =
      !execution.coefficientRestartFamilies.empty() &&
      execution.outputFields == observer.stateBlockIdentifiers &&
      execution.coefficientRestartFamilies ==
          observer.stateBlockIdentifiers &&
      observer.fieldNames.empty() && observer.x.empty() && observer.y.empty() &&
      observer.z.empty();
  if (isCompleteCoefficientObserver) {
    const auto *block = observer.stateBlockIdentifiers.empty()
                            ? nullptr
                            : context.stateBlock(
                                  observer.stateBlockIdentifiers.front());
    plan.schema.metadata.variables.push_back(
        metadataReal("absTolerance",
                     block == nullptr ? 1e-6 : block->absoluteTolerance));
  }

  WVLegacyObserverOperationBinder operations;
  operations.fullField = [&]() -> WVKernelStatus {
    for (const auto &field : execution.outputFields) {
      const auto status = addFullField(field);
      if (!status)
        return status;
    }
    return WVKernelStatus::ok();
  };
  operations.fixedVerticalProfiles = [&]() -> WVKernelStatus {
    if (configuration == nullptr || spatialDimensions.size() != 3)
      return {WVKernelStatusCode::unsupportedOperation,
              "Fixed vertical profiles require a three-dimensional transform."};
    const auto &legacyConfiguration = *configuration;
    const std::string idName = observer.name + "_id";
    const std::string zName = observer.name + "_z";
    addAxis(plan.schema, idName, observer.x.size(),
            WVObservationCoordinateRole::identifier);
    addAxis(plan.schema, zName, legacyConfiguration.Nz,
            WVObservationCoordinateRole::z);
    std::vector<double> identifiers(observer.x.size());
    for (std::size_t index = 0; index < identifiers.size(); ++index)
      identifiers[index] = static_cast<double>(index + 1);
    addConstantReal(plan, "static-" + idName, idName, {idName},
                    WVObservationValueLayout::staticValue,
                    std::move(identifiers), "unitless id number", "",
                    WVObservationCoordinateRole::identifier);
    std::vector<double> z = observer.z;
    if (z.empty()) {
      z.resize(legacyConfiguration.Nz);
      const double dz =
          legacyConfiguration.Lz /
          static_cast<double>(legacyConfiguration.Nz - 1);
      for (std::size_t index = 0; index < z.size(); ++index)
        z[index] = -legacyConfiguration.Lz +
                   static_cast<double>(index) * dz;
    }
    addConstantReal(plan, "static-" + zName, zName, {zName},
                    WVObservationValueLayout::staticValue, std::move(z), "m",
                    "z-positions of mooring observations",
                    WVObservationCoordinateRole::z);
    std::vector<double> x = observer.x;
    std::vector<double> y = observer.y;
    WVFieldSamplingRequest sampling;
    sampling.kind = WVFieldSamplingKind::fixedVerticalProfiles;
    const double dx = legacyConfiguration.Lx /
                      static_cast<double>(legacyConfiguration.Nx);
    const double dy = legacyConfiguration.Ly /
                      static_cast<double>(legacyConfiguration.Ny);
    for (std::size_t index = 0; index < x.size(); ++index) {
      x[index] = std::fmod(x[index], legacyConfiguration.Lx);
      y[index] = std::fmod(y[index], legacyConfiguration.Ly);
      if (x[index] < 0.0)
        x[index] += legacyConfiguration.Lx;
      if (y[index] < 0.0)
        y[index] += legacyConfiguration.Ly;
      sampling.xIndices.push_back(std::min(
          legacyConfiguration.Nx,
          static_cast<std::size_t>(std::floor(x[index] / dx)) + 1));
      sampling.yIndices.push_back(std::min(
          legacyConfiguration.Ny,
          static_cast<std::size_t>(std::floor(y[index] / dy)) + 1));
    }
    addConstantReal(plan, "static-x", observer.name + "_x", {idName},
                    WVObservationValueLayout::staticValue, std::move(x), "m",
                    "x position of mooring", WVObservationCoordinateRole::x);
    addConstantReal(plan, "static-y", observer.name + "_y", {idName},
                    WVObservationValueLayout::staticValue, std::move(y), "m",
                    "y position of mooring", WVObservationCoordinateRole::y);
    for (const auto &field : execution.outputFields) {
      const auto *metadata = findPortableVariable(field);
      if (metadata == nullptr ||
          metadata->kind != WVPortableVariableKind::field)
        return invalid("Unsupported mooring field: " + field + ".");
      WVObserverOutputChannel channel;
      channel.source = WVObserverOutputChannelSource::sampledField;
      channel.sourceIdentifier = field;
      channel.sampling = sampling;
      addChannel(
          plan,
          fieldVariable(*metadata, "derived-" + field,
                        observer.name + '_' + field, {zName, idName},
                        outputLayout(*metadata), ", recorded at the mooring"),
          std::move(channel));
    }
    return WVKernelStatus::ok();
  };
  operations.fixedPositions = [&]() -> WVKernelStatus {
    plan.schema.metadata.variables.push_back(
        metadataReal("outputScale", observer.outputScale));
    plan.schema.metadata.variables.push_back(
        metadataReal("outputOffset", observer.outputOffset));
    plan.schema.metadata.attributes.push_back(
        {"trackedVarInterpolation",
         observer.trackedFieldInterpolation == WVPositionInterpolation::linear
             ? "linear"
             : "spline"});
    const std::string idName = observer.name + "_id";
    addAxis(plan.schema, idName, observer.x.size(),
            WVObservationCoordinateRole::identifier);
    std::vector<double> identifiers(observer.x.size());
    for (std::size_t index = 0; index < identifiers.size(); ++index)
      identifiers[index] = static_cast<double>(index + 1);
    addConstantReal(plan, "static-" + idName, idName, {idName},
                    WVObservationValueLayout::staticValue,
                    std::move(identifiers), "unitless id number", "",
                    WVObservationCoordinateRole::identifier);
    for (const auto &[suffix, values, role] :
         std::array<std::tuple<const char *, const std::vector<double> *,
                               WVObservationCoordinateRole>,
                    3>{{{"x", &observer.x, WVObservationCoordinateRole::x},
                        {"y", &observer.y, WVObservationCoordinateRole::y},
                        {"z", &observer.z, WVObservationCoordinateRole::z}}})
      addConstantReal(plan, "static-" + std::string(suffix),
                      observer.name + '_' + suffix, {idName},
                      WVObservationValueLayout::staticValue, *values, "m",
                      std::string(suffix) + " position of fixed observation",
                      role);
    const auto &field = execution.outputFields.front();
    const auto *metadata = findPortableVariable(field);
    if (metadata == nullptr ||
        metadata->kind != WVPortableVariableKind::field)
      return invalid("Unsupported fixed-position field: " + field + ".");
    WVObserverOutputChannel channel;
    channel.source = WVObserverOutputChannelSource::sampledField;
    channel.sourceIdentifier = field;
    channel.sampling.kind = WVFieldSamplingKind::positions;
    channel.sampling.x = observer.x;
    channel.sampling.y = observer.y;
    channel.sampling.z = observer.z;
    channel.sampling.interpolation = observer.trackedFieldInterpolation;
    channel.scale = observer.outputScale;
    channel.offset = observer.outputOffset;
    addChannel(
        plan,
        fieldVariable(*metadata, "derived-" + field,
                      observer.name + "_value", {idName},
                      outputLayout(*metadata),
                      ", sampled and affinely transformed by the observing system"),
        std::move(channel));
    return WVKernelStatus::ok();
  };
  operations.movingPositions = [&]() -> WVKernelStatus {
    plan.schema.metadata.variables.push_back(
        metadataBoolean("isXYOnly", observer.isXYOnly));
    plan.schema.metadata.variables.push_back(metadataReal(
        "absToleranceXY", observer.horizontalAbsoluteTolerance));
    plan.schema.metadata.variables.push_back(metadataReal(
        "absToleranceZ", observer.verticalAbsoluteTolerance));
    plan.schema.metadata.attributes.push_back(
        {"advectionInterpolation",
         observer.advectionInterpolation == WVPositionInterpolation::linear
             ? "linear"
             : "spline"});
    plan.schema.metadata.attributes.push_back(
        {"trackedVarInterpolation",
         observer.trackedFieldInterpolation == WVPositionInterpolation::linear
             ? "linear"
             : "spline"});
    const std::string idName = observer.name + "_id";
    addAxis(plan.schema, idName, observer.x.size(),
            WVObservationCoordinateRole::identifier);
    std::vector<double> identifiers(observer.x.size());
    for (std::size_t index = 0; index < identifiers.size(); ++index)
      identifiers[index] = static_cast<double>(index + 1);
    const auto particleAttributes = [&](const std::string &value) {
      return std::vector<WVObservationAttribute>{
          {"isParticle", "1"}, {"particleName", observer.name},
          {"particleVariableName", value}};
    };
    addConstantReal(plan, "static-" + idName, idName, {idName},
                    WVObservationValueLayout::staticValue,
                    std::move(identifiers), "unitless id number", "",
                    WVObservationCoordinateRole::identifier,
                    particleAttributes("id"));
    const auto channels = particlePositionChannels(observer.isXYOnly);
    for (std::size_t index = 0; index < channels.size(); ++index) {
      const std::string suffix = movingFieldChannelName(channels[index]);
      WVObservationVariable variable;
      variable.identifier = "state-" + suffix;
      variable.name = movingFieldVariableName(observer, channels[index]);
      variable.dimensionIdentifiers = {idName};
      variable.layout = WVObservationValueLayout::record;
      variable.units = "m";
      variable.description = suffix + " position of particle";
      variable.attributes = particleAttributes(suffix);
      variable.coordinateRole =
          channels[index] == WVMovingFieldChannel::x
              ? WVObservationCoordinateRole::x
              : channels[index] == WVMovingFieldChannel::y
                    ? WVObservationCoordinateRole::y
                    : WVObservationCoordinateRole::z;
      WVObserverOutputChannel channel;
      channel.source = WVObserverOutputChannelSource::additionalState;
      channel.sourceIdentifier = observer.stateBlockIdentifiers[index];
      addChannel(plan, std::move(variable), std::move(channel));
    }
    if (observer.isXYOnly && !observer.z.empty())
      addConstantReal(plan, "fixed-z", observer.name + "_z", {idName},
                      WVObservationValueLayout::record, observer.z, "m",
                      "z position of particle",
                      WVObservationCoordinateRole::z,
                      particleAttributes("z"));
    plan.movingPositions.stateBlockIdentifiers =
        observer.stateBlockIdentifiers;
    plan.movingPositions.fixedZ = observer.z;
    plan.movingPositions.positionCount = observer.x.size();
    plan.movingPositions.isXYOnly = observer.isXYOnly;
    plan.movingPositions.interpolation =
        observer.trackedFieldInterpolation;
    for (const auto &field : execution.outputFields) {
      const auto *metadata = findPortableVariable(field);
      if (metadata == nullptr ||
          metadata->kind != WVPortableVariableKind::field)
        return invalid("Unsupported particle tracked field: " + field + ".");
      auto variable = fieldVariable(
          *metadata, "derived-" + field, observer.name + '_' + field,
          {idName}, WVObservationValueLayout::record,
          ", recorded along the particle trajectory");
      variable.attributes = particleAttributes(field);
      WVObserverOutputChannel channel;
      channel.source = WVObserverOutputChannelSource::movingField;
      channel.sourceIdentifier = field;
      addChannel(plan, std::move(variable), std::move(channel));
    }
    return WVKernelStatus::ok();
  };
  operations.integratedState = [&]() -> WVKernelStatus {
    plan.schema.metadata.variables.push_back(
        metadataBoolean("isXYOnly", observer.isXYOnly));
    const auto *block =
        context.stateBlock(observer.stateBlockIdentifiers.front());
    if (block == nullptr)
      return invalid("Integrated observer state block is unavailable.");
    plan.schema.metadata.variables.push_back(metadataReal(
        "absTolerance", block->absoluteTolerance));
    plan.schema.metadata.variables.push_back(
        metadataBoolean("shouldAntialias", observer.shouldAntialias));
    const std::vector<std::string> names =
        block->dimensions.size() == 2
            ? std::vector<std::string>{"x", "y"}
            : std::vector<std::string>{"x", "y", "z"};
    for (std::size_t index = 0; index < names.size(); ++index)
      addAxis(plan.schema, names[index], block->dimensions[index],
              names[index] == "x"
                  ? WVObservationCoordinateRole::x
                  : names[index] == "y" ? WVObservationCoordinateRole::y
                                           : WVObservationCoordinateRole::z);
    WVObservationVariable variable;
    variable.identifier = "state-tracer";
    variable.name = observer.name;
    variable.dimensionIdentifiers = names;
    variable.layout = WVObservationValueLayout::record;
    variable.attributes = {{"isTracer", "1"}};
    WVObserverOutputChannel channel;
    channel.source = WVObserverOutputChannelSource::additionalState;
    channel.sourceIdentifier = observer.stateBlockIdentifiers.front();
    addChannel(plan, std::move(variable), std::move(channel));
    return WVKernelStatus::ok();
  };
  const auto operationStatus = resolveOperation(observer, operations);
  if (!operationStatus)
    return operationStatus;
  const auto status = validateObservationSchema(plan.schema);
  if (!status)
    return status;
  output = std::move(plan);
  return WVKernelStatus::ok();
}

using BuiltInExecutionPlanResolver = WVKernelStatus (*)(
    const WVObserverRecord &, WVObserverExecutionPlan &);

WVKernelStatus coefficientExecutionPlan(const WVObserverRecord &observer,
                                        WVObserverExecutionPlan &plan) {
  plan = {};
  plan.persistedName = observer.name;
  plan.outputFields = observer.stateBlockIdentifiers;
  plan.coefficientRestartFamilies = plan.outputFields;
  return WVKernelStatus::ok();
}

WVKernelStatus eulerianExecutionPlan(const WVObserverRecord &observer,
                                     WVObserverExecutionPlan &plan) {
  plan = {};
  plan.fieldListAttribute = "fieldNames";
  auto status = configuredText(observer, "fieldNames", observer.fieldNames,
                               plan.outputFields);
  if (!status)
    return status;
  for (const auto *family : {"Ap", "Am", "A0"})
    if (std::find(observer.fieldNames.begin(), observer.fieldNames.end(),
                  family) != observer.fieldNames.end())
      plan.coefficientRestartFamilies.emplace_back(family);
  return WVKernelStatus::ok();
}

WVKernelStatus trackedFieldExecutionPlan(const WVObserverRecord &observer,
                                         WVObserverExecutionPlan &plan) {
  plan = {};
  plan.fieldListAttribute = "trackedFieldNames";
  plan.persistedName = observer.name;
  return configuredText(observer, "trackedFieldNames", observer.fieldNames,
                        plan.outputFields);
}

WVKernelStatus tracerExecutionPlan(const WVObserverRecord &observer,
                                   WVObserverExecutionPlan &plan) {
  plan = {};
  plan.persistedName = observer.name;
  return WVKernelStatus::ok();
}

class BuiltInObservingSystem : public WVObservingSystem {
public:
  BuiltInObservingSystem(std::string typeIdentifier,
                         WVPortableTypedRecord configuration,
                         WVLegacyObserverOperationResolver resolveOperation)
      : typeIdentifier_(std::move(typeIdentifier)),
        configuration_(std::move(configuration)),
        resolveOperation_(std::move(resolveOperation)) {}

  const std::string &typeIdentifier() const noexcept override {
    return typeIdentifier_;
  }
  std::uint32_t contractVersion() const noexcept override {
    return WVPortablePairContractVersion;
  }
  WVKernelStatus outputPlan(
      const WVObserverRecord &observer,
      const WVObserverOutputPlanningContext &context,
      WVObserverOutputPlan &plan) const override {
    WVObserverExecutionPlan execution;
    const auto status = executionPlan(observer, execution);
    if (!status)
      return status;
    return buildLegacyOutputPlan(typeIdentifier(), observer, execution,
                                 resolveOperation_, context, plan);
  }
  std::size_t persistentBytes() const noexcept override {
    return sizeof(*this) + typeIdentifier_.capacity() +
           configuration_.persistentBytes() - sizeof(configuration_);
  }

protected:
  static WVKernelStatus validateSampleOnly(
      const WVObserverRecord &observer) {
    return observer.stateBlockIdentifiers.empty()
               ? WVKernelStatus::ok()
               : invalid("Sample-only observers cannot own integrated state blocks.");
  }

private:
  std::string typeIdentifier_;
  WVPortableTypedRecord configuration_;
  WVLegacyObserverOperationResolver resolveOperation_;
};

class WVCoefficientsImplementation final : public BuiltInObservingSystem {
public:
  explicit WVCoefficientsImplementation(WVPortableTypedRecord configuration)
      : BuiltInObservingSystem("WVCoefficients", std::move(configuration),
                               &selectFullFieldOperation) {}
  WVKernelStatus executionPlan(const WVObserverRecord &observer,
                               WVObserverExecutionPlan &plan) const override {
    return coefficientExecutionPlan(observer, plan);
  }
  WVKernelStatus validate(
      const WVObserverRecord &observer,
      const std::map<std::string, const WVStateBlockRecord *> &,
      std::map<std::string, std::size_t> &) const override {
    if (observer.stateBlockIdentifiers.empty())
      return invalid("WVCoefficients must reference at least one coefficient "
                     "family in transform order.");
    return WVKernelStatus::ok();
  }
};

class WVEulerianFieldsImplementation final : public BuiltInObservingSystem {
public:
  explicit WVEulerianFieldsImplementation(WVPortableTypedRecord configuration)
      : BuiltInObservingSystem("WVEulerianFields", std::move(configuration),
                               &selectFullFieldOperation) {}
  WVKernelStatus executionPlan(const WVObserverRecord &observer,
                               WVObserverExecutionPlan &plan) const override {
    return eulerianExecutionPlan(observer, plan);
  }
  WVKernelStatus validate(
      const WVObserverRecord &observer,
      const std::map<std::string, const WVStateBlockRecord *> &,
      std::map<std::string, std::size_t> &) const override {
    return validateSampleOnly(observer);
  }
};

class WVMooringImplementation final : public BuiltInObservingSystem {
public:
  explicit WVMooringImplementation(WVPortableTypedRecord configuration)
      : BuiltInObservingSystem("WVMooring", std::move(configuration),
                               &selectFixedVerticalProfileOperation) {}
  WVKernelStatus executionPlan(const WVObserverRecord &observer,
                               WVObserverExecutionPlan &plan) const override {
    return trackedFieldExecutionPlan(observer, plan);
  }
  WVKernelStatus validate(
      const WVObserverRecord &observer,
      const std::map<std::string, const WVStateBlockRecord *> &,
      std::map<std::string, std::size_t> &) const override {
    const auto status = validateSampleOnly(observer);
    if (!status)
      return status;
    if (observer.x.empty() || observer.x.size() != observer.y.size())
      return invalid("WVMooring requires equal nonempty x and y coordinates.");
    return WVKernelStatus::ok();
  }
};

class WVLagrangianParticlesImplementation final
    : public BuiltInObservingSystem {
public:
  explicit WVLagrangianParticlesImplementation(
      WVPortableTypedRecord configuration)
      : BuiltInObservingSystem("WVLagrangianParticles",
                               std::move(configuration),
                               &selectMovingPositionOperation) {}
  WVKernelStatus executionPlan(const WVObserverRecord &observer,
                               WVObserverExecutionPlan &plan) const override {
    return trackedFieldExecutionPlan(observer, plan);
  }
  WVKernelStatus bindIntegration(
      const WVObserverRecord &observer,
      const WVObserverIntegrationBinder &binder) const override {
    if (!binder.advectedPositions)
      return invalid("The integration runtime cannot bind advected positions.");
    return binder.advectedPositions(observer);
  }
  WVKernelStatus validate(
      const WVObserverRecord &observer,
      const std::map<std::string, const WVStateBlockRecord *> &blocks,
      std::map<std::string, std::size_t> &ownerCounts) const override {
    if (observer.x.empty() || observer.x.size() != observer.y.size())
      return invalid(
          "WVLagrangianParticles requires equal nonempty x and y coordinates.");
    if (!observer.isXYOnly && observer.z.size() != observer.x.size())
      return invalid(
          "Three-dimensional particles require one z coordinate per particle.");
    if (!(observer.horizontalAbsoluteTolerance > 0.0) ||
        !std::isfinite(observer.horizontalAbsoluteTolerance))
      return invalid(
          "Particle horizontal tolerance must be finite and positive.");
    if (!observer.isXYOnly &&
        (!(observer.verticalAbsoluteTolerance > 0.0) ||
         !std::isfinite(observer.verticalAbsoluteTolerance)))
      return invalid("Particle vertical tolerance must be finite and positive.");
    const auto channels = particlePositionChannels(observer.isXYOnly);
    if (observer.stateBlockIdentifiers.size() != channels.size())
      return invalid("WVLagrangianParticles requires ordered x, y, and "
                     "optional z state blocks.");
    for (const auto &identifier : observer.stateBlockIdentifiers) {
      const auto found = blocks.find(identifier);
      if (found == blocks.end())
        return invalid("Particle observer references an unknown state block.");
      const auto *block = found->second;
      if (block->scalarType != WVStateScalarType::real64 ||
          block->ownership != WVStateOwnership::integratorOwned ||
          block->dimensions !=
              std::vector<std::size_t>({observer.x.size()}))
        return invalid("Particle state blocks must be integrator-owned real "
                       "vectors matching the particle count.");
      ++ownerCounts.at(identifier);
    }
    return WVKernelStatus::ok();
  }
};

class WVTracerImplementation final : public BuiltInObservingSystem {
public:
  explicit WVTracerImplementation(WVPortableTypedRecord configuration)
      : BuiltInObservingSystem("WVTracer", std::move(configuration),
                               &selectIntegratedStateOperation) {}
  WVKernelStatus executionPlan(const WVObserverRecord &observer,
                               WVObserverExecutionPlan &plan) const override {
    return tracerExecutionPlan(observer, plan);
  }
  WVKernelStatus bindIntegration(
      const WVObserverRecord &observer,
      const WVObserverIntegrationBinder &binder) const override {
    if (!binder.advectedScalar)
      return invalid("The integration runtime cannot bind an advected scalar.");
    return binder.advectedScalar(observer);
  }
  WVKernelStatus validate(
      const WVObserverRecord &observer,
      const std::map<std::string, const WVStateBlockRecord *> &blocks,
      std::map<std::string, std::size_t> &ownerCounts) const override {
    if (observer.stateBlockIdentifiers.size() != 1)
      return invalid("WVTracer requires exactly one state block.");
    const auto found = blocks.find(observer.stateBlockIdentifiers.front());
    if (found == blocks.end())
      return invalid("WVTracer references an unknown state block.");
    const auto *block = found->second;
    if (block->scalarType != WVStateScalarType::real64 ||
        block->ownership != WVStateOwnership::integratorOwned)
      return invalid("WVTracer requires one integrator-owned real state block.");
    const std::size_t expectedRank = observer.isXYOnly ? 2 : 3;
    if (block->dimensions.size() != expectedRank)
      return invalid(observer.isXYOnly
                         ? "A two-dimensional WVTracer requires a rank-two state block."
                         : "A three-dimensional WVTracer requires a rank-three state block.");
    ++ownerCounts.at(observer.stateBlockIdentifiers.front());
    return WVKernelStatus::ok();
  }
};

} // namespace

WVKernelStatus addBuiltInObserverFactories(
    WVExtensionCatalogBuilder &builder) {
  const auto add = [&](std::string identity, auto constructor,
                       LegacyConfigurationPopulator populate,
                       WVLegacyObserverOperationResolver resolveOperation,
                       BuiltInExecutionPlanResolver resolveExecutionPlan,
                       WVLegacyObserverPersistenceMetadata persistence) {
    WVObserverOutputPlanResolver resolveOutputPlan =
        [identity, resolveOperation, resolveExecutionPlan](
            const WVObserverRecord &observer,
            const WVObserverOutputPlanningContext &context,
            WVObserverOutputPlan &plan) {
          WVObserverExecutionPlan execution;
          auto status = resolveExecutionPlan(observer, execution);
          if (!status)
            return status;
          return buildLegacyOutputPlan(identity, observer, execution,
                                       resolveOperation, context, plan);
        };
    WVObserverFactoryRegistration registration(
        std::move(identity), WVPortablePairContractVersion,
        [constructor](const WVObserverRecord &,
                      const WVPortableTypedRecord &configuration,
                      std::shared_ptr<const WVObservingSystem> &result) {
          try {
            result = constructor(configuration);
            return WVKernelStatus::ok();
          } catch (const std::bad_alloc &) {
            return WVKernelStatus{WVKernelStatusCode::allocationFailure,
                                  "Unable to construct an observer implementation."};
          }
        },
        [populate](const WVObserverRecord &observer,
                   WVPortableTypedRecord &configuration) {
          return resolveLegacyConfiguration(observer, populate,
                                            configuration);
        },
        std::move(resolveOutputPlan));
    auto status =
        WVObserverFactoryRegistrationAccess::attachLegacyCompatibility(
            registration, std::move(resolveOperation),
            std::move(persistence));
    if (!status)
      return status;
    return builder.addObserverFactory(std::move(registration));
  };
  auto status = add("WVCoefficients", [](const auto &configuration) {
    return std::make_shared<WVCoefficientsImplementation>(configuration);
  }, &populateCoefficientConfiguration, &selectFullFieldOperation,
  &coefficientExecutionPlan,
  {{}, {"Ap", "Am", "A0"}, "coefficients", false});
  if (!status) return status;
  status = add("WVEulerianFields", [](const auto &configuration) {
    return std::make_shared<WVEulerianFieldsImplementation>(configuration);
  }, &populateEulerianConfiguration, &selectFullFieldOperation,
  &eulerianExecutionPlan,
  {"fieldNames", {}, "eulerian-fields", true});
  if (!status) return status;
  status = add("WVMooring", [](const auto &configuration) {
    return std::make_shared<WVMooringImplementation>(configuration);
  }, &populateMooringConfiguration, &selectFixedVerticalProfileOperation,
  &trackedFieldExecutionPlan,
  {"trackedFieldNames", {}, {}, false});
  if (!status) return status;
  status = add("WVLagrangianParticles", [](const auto &configuration) {
    return std::make_shared<WVLagrangianParticlesImplementation>(configuration);
  }, &populateParticleConfiguration, &selectMovingPositionOperation,
  &trackedFieldExecutionPlan,
  {"trackedFieldNames", {}, {}, false});
  if (!status) return status;
  return add("WVTracer", [](const auto &configuration) {
    return std::make_shared<WVTracerImplementation>(configuration);
  }, &populateTracerConfiguration, &selectIntegratedStateOperation,
  &tracerExecutionPlan,
  {});
}

WVKernelStatus canonicalCoefficientObserver(
    std::string identifier, const WVExtensionCatalog &catalog,
    WVObserverRecord &observer,
    std::vector<std::string> coefficientFamilies) {
  if (catalog.observers().registration(
          "WVCoefficients", WVPortablePairContractVersion) == nullptr)
    return {WVKernelStatusCode::unsupportedOperation,
            "No coefficient observer implementation is registered."};
  WVObserverRecord candidate;
  candidate.identifier = std::move(identifier);
  candidate.name = "WVCoefficients";
  candidate.typeIdentifier = "WVCoefficients";
  candidate.contractVersion = WVPortablePairContractVersion;
  if (coefficientFamilies.empty())
    return invalid("A canonical coefficient observer requires at least one "
                   "resolved transform family.");
  candidate.stateBlockIdentifiers = std::move(coefficientFamilies);
  const auto configurationStatus =
      catalog.observers().resolveConfiguration(candidate,
                                               candidate.configuration);
  if (!configurationStatus)
    return configurationStatus;
  observer = std::move(candidate);
  return WVKernelStatus::ok();
}

const char *movingFieldChannelName(WVMovingFieldChannel channel) noexcept {
  switch (channel) {
  case WVMovingFieldChannel::x:
    return "x";
  case WVMovingFieldChannel::y:
    return "y";
  case WVMovingFieldChannel::z:
    return "z";
  case WVMovingFieldChannel::tracerValue:
    return "value";
  }
  return nullptr;
}

std::vector<WVMovingFieldChannel> particlePositionChannels(bool isXYOnly) {
  return isXYOnly
             ? std::vector<WVMovingFieldChannel>{WVMovingFieldChannel::x,
                                                 WVMovingFieldChannel::y}
             : std::vector<WVMovingFieldChannel>{WVMovingFieldChannel::x,
                                                 WVMovingFieldChannel::y,
                                                 WVMovingFieldChannel::z};
}

std::string movingFieldVariableName(const WVObserverRecord &observer,
                                    WVMovingFieldChannel channel) {
  if (channel == WVMovingFieldChannel::tracerValue)
    return observer.name;
  return observer.name + '_' + movingFieldChannelName(channel);
}

} // namespace wavevortex::runtime::detail
