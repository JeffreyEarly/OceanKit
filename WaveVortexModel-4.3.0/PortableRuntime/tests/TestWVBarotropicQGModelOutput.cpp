#include "WaveVortexRuntime/WVExtensionCatalog.hpp"
#include "WaveVortexRuntime/WVModel.hpp"
#include "WVReferenceFFTEngine.hpp"

#include <netcdf.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace wavevortex;
using namespace wavevortex::runtime;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

struct TemporaryDirectory {
  std::filesystem::path path =
      std::filesystem::temp_directory_path() /
      ("wave-vortex-qg-model-output-" +
       std::to_string(
           std::chrono::steady_clock::now().time_since_epoch().count()));
  TemporaryDirectory() { std::filesystem::create_directories(path); }
  ~TemporaryDirectory() {
    std::error_code ignored;
    std::filesystem::remove_all(path, ignored);
  }
};

std::shared_ptr<const WVExtensionCatalog> catalog() {
  std::shared_ptr<const WVExtensionCatalog> result;
  const auto status = makeBuiltInExtensionCatalog(result);
  require(static_cast<bool>(status), status.message);
  return result;
}

WVTransformBarotropicQGConfiguration qgConfiguration() {
  WVTransformBarotropicQGConfiguration configuration;
  configuration.Nx = 5;
  configuration.Ny = 4;
  configuration.Lx = 15000.0;
  configuration.Ly = 9000.0;
  configuration.h = 0.8;
  configuration.j = 1;
  configuration.g = 9.80665;
  configuration.planetaryRadius = 6.3712e6;
  configuration.rotationRate = 7.292115e-5;
  configuration.latitude = 33.0;
  configuration.shouldAntialias = true;
  return configuration;
}

WVTransformStateDescription descriptionFor(
    const WVIntegrationStateLayout &layout) {
  WVTransformStateDescription result;
  result.transformIdentifier = layout.transformIdentifier();
  result.spatialDimensions = layout.spatialDimensions();
  result.supportsFixedTimeStepSelection = true;
  for (const auto &family : layout.coefficientFamilies())
    result.coefficientFamilies.push_back(
        {family.identifier, family.spectralDimensions, family.toleranceKind});
  return result;
}

WVPortableObserverRecord observerRecord(
    const WVTransformBarotropicQGConfiguration &configuration,
    std::size_t Nkl, const std::filesystem::path &path, double finalTime) {
  WVPortableObserverRecord record;
  record.stateBlocks = {
      {"A0", WVStateScalarType::complex64, {Nkl},
       WVToleranceKind::coefficientEnergyScaled, 1e-8,
       WVStateOwnership::integratorOwned,
       WVRestartRequirement::requiredDynamicState},
      {"floats-x", WVStateScalarType::real64, {2},
       WVToleranceKind::uniformAbsolute, 1e-5,
       WVStateOwnership::integratorOwned,
       WVRestartRequirement::requiredDynamicState},
      {"floats-y", WVStateScalarType::real64, {2},
       WVToleranceKind::uniformAbsolute, 1e-5,
       WVStateOwnership::integratorOwned,
       WVRestartRequirement::requiredDynamicState},
      {"dye-state", WVStateScalarType::real64,
       {configuration.Nx, configuration.Ny},
       WVToleranceKind::uniformAbsolute, 1e-7,
       WVStateOwnership::integratorOwned,
       WVRestartRequirement::requiredDynamicState}};

  WVObserverRecord coefficients;
  coefficients.identifier = "coefficients";
  coefficients.name = "WVCoefficients";
  coefficients.typeIdentifier = "WVCoefficients";
  coefficients.stateBlockIdentifiers = {"A0"};

  WVObserverRecord fields;
  fields.identifier = "fields";
  fields.name = "WVEulerianFields";
  fields.typeIdentifier = "WVEulerianFields";
  fields.fieldNames = {"A0", "u", "v", "eta", "pi", "psi", "qgpv",
                       "zeta_z", "ssh", "energy", "uvMax"};

  WVObserverRecord particles;
  particles.identifier = "floats";
  particles.name = "floats";
  particles.typeIdentifier = "WVLagrangianParticles";
  particles.stateBlockIdentifiers = {"floats-x", "floats-y"};
  particles.fieldNames = {"u", "qgpv"};
  particles.x = {0.13 * configuration.Lx, 0.77 * configuration.Lx};
  particles.y = {0.21 * configuration.Ly, 0.63 * configuration.Ly};
  particles.isXYOnly = true;
  particles.horizontalAbsoluteTolerance = 1e-5;
  particles.advectionInterpolation = WVPositionInterpolation::linear;
  particles.trackedFieldInterpolation = WVPositionInterpolation::spline;

  WVObserverRecord tracer;
  tracer.identifier = "dye";
  tracer.name = "dye";
  tracer.typeIdentifier = "WVTracer";
  tracer.stateBlockIdentifiers = {"dye-state"};
  tracer.isXYOnly = true;
  tracer.shouldAntialias = true;

  record.observers = {coefficients, fields, particles, tracer};
  record.outputFiles = {
      {"primary",
       path.string(),
       {{"restart", "qg",
         {0.01, 0.0, finalTime},
         {"coefficients", "fields", "floats", "dye"}, true}}}};
  return record;
}

WVCheckpoint checkpointFor(
    const WVTransformBarotropicQGConfiguration &configuration,
    const WVTransformStateDescription &description,
    const WVFrozenForcingSchedule &forcingSchedule) {
  WVCheckpoint checkpoint;
  checkpoint.transformKind = WVPersistedTransformKind::barotropicQG;
  checkpoint.barotropicQGConfiguration = configuration;
  checkpoint.stateDescription = description;
  checkpoint.state.t = 0.0;
  checkpoint.state.t0 = 0.0;
  checkpoint.transformState.transformIdentifier =
      description.transformIdentifier;
  checkpoint.transformState.spatialDimensions =
      description.spatialDimensions;
  checkpoint.transformState.t = 0.0;
  checkpoint.transformState.t0 = 0.0;
  for (const auto &family : description.coefficientFamilies) {
    WVCoefficientFamilyCheckpoint values;
    values.identifier = family.identifier;
    values.spectralDimensions = family.spectralDimensions;
    std::size_t count = 1;
    for (const auto extent : family.spectralDimensions)
      count *= extent;
    values.values.resize(count);
    for (std::size_t index = 0; index < count; ++index) {
      const double p = static_cast<double>(index + 1);
      values.values[index] = {2e-5 * std::sin(0.37 * p),
                              1e-5 * std::cos(0.23 * (p + 2.0))};
    }
    checkpoint.transformState.coefficientFamilies.push_back(
        std::move(values));
  }
  checkpoint.metadata.modelVersion = "4.3.0";
  checkpoint.metadata.transformClass = "WVTransformBarotropicQG";
  checkpoint.metadata.stateGroupPath = "/qg";
  checkpoint.metadata.stateCount = 1;
  checkpoint.forcingSchedule = forcingSchedule;
  return checkpoint;
}

void initializeAdditionalState(WVModel &model, WVModelState &state) {
  auto status = model.initializeObserverState(state);
  require(static_cast<bool>(status), status.message);
  auto view = state.mutableView();
  for (std::size_t block = 0; block < view.additionalBlockCount; ++block) {
    auto &values = view.additionalBlocks[block];
    if (values.layout->identifier != "dye-state")
      continue;
    const auto &dimensions = values.layout->dimensions;
    for (std::size_t iy = 0; iy < dimensions[1]; ++iy)
      for (std::size_t ix = 0; ix < dimensions[0]; ++ix)
        values.realData[ix + dimensions[0] * iy] =
            0.3 + 0.2 * std::sin(0.7 * static_cast<double>(ix)) *
                      std::cos(0.4 * static_cast<double>(iy));
  }
  status = model.prepareStateAfterRestart(state);
  require(static_cast<bool>(status), status.message);
}

double relativeStateDifference(WVModelState &left, WVModelState &right) {
  const auto a = left.constView();
  const auto b = right.constView();
  require(a.coefficientFamilyCount == b.coefficientFamilyCount &&
              a.additionalBlockCount == b.additionalBlockCount,
          "state layouts differ");
  double maximum = 0.0;
  for (std::size_t family = 0; family < a.coefficientFamilyCount; ++family) {
    require(a.coefficientFamilies[family].layout->identifier ==
                b.coefficientFamilies[family].layout->identifier &&
                a.coefficientFamilies[family].layout->elementCount ==
                    b.coefficientFamilies[family].layout->elementCount,
            "coefficient families differ");
    for (std::size_t index = 0;
         index < a.coefficientFamilies[family].layout->elementCount; ++index) {
      const auto av = a.coefficientFamilies[family].data[index];
      const auto bv = b.coefficientFamilies[family].data[index];
      maximum = std::max(
          maximum,
          std::hypot(av.real - bv.real, av.imag - bv.imag) /
              std::max(1.0, std::hypot(av.real, av.imag)));
    }
  }
  for (std::size_t block = 0; block < a.additionalBlockCount; ++block) {
    require(a.additionalBlocks[block].layout->identifier ==
                b.additionalBlocks[block].layout->identifier &&
                a.additionalBlocks[block].layout->elementCount ==
                    b.additionalBlocks[block].layout->elementCount,
            "additional state blocks differ");
    for (std::size_t index = 0;
         index < a.additionalBlocks[block].layout->elementCount; ++index) {
      const double av = a.additionalBlocks[block].realData[index];
      const double bv = b.additionalBlocks[block].realData[index];
      maximum = std::max(maximum,
                         std::abs(av - bv) / std::max(1.0, std::abs(av)));
    }
  }
  return maximum;
}

bool findVariable(int group, const std::string &name, int &owner, int &variable) {
  if (nc_inq_varid(group, name.c_str(), &variable) == NC_NOERR) {
    owner = group;
    return true;
  }
  int count = 0;
  if (nc_inq_grps(group, &count, nullptr) != NC_NOERR || count == 0)
    return false;
  std::vector<int> groups(static_cast<std::size_t>(count));
  require(nc_inq_grps(group, &count, groups.data()) == NC_NOERR,
          "NetCDF group lookup failed");
  for (const auto child : groups)
    if (findVariable(child, name, owner, variable))
      return true;
  return false;
}

void verifyCompactSchema(const std::filesystem::path &path) {
  int file = -1;
  require(nc_open(path.c_str(), NC_NOWRITE, &file) == NC_NOERR,
          "unable to inspect QG output schema");
  int owner = -1;
  int variable = -1;
  require(findVariable(file, "A0_real", owner, variable),
          "compact A0_real output is absent");
  int rank = 0;
  require(nc_inq_varndims(owner, variable, &rank) == NC_NOERR && rank == 2,
          "compact A0 output must have [t,kl] rank");
  int dimensions[NC_MAX_VAR_DIMS] = {};
  require(nc_inq_vardimid(owner, variable, dimensions) == NC_NOERR,
          "unable to inspect compact A0 dimensions");
  char firstName[NC_MAX_NAME + 1] = {};
  char secondName[NC_MAX_NAME + 1] = {};
  require(nc_inq_dimname(owner, dimensions[0], firstName) == NC_NOERR &&
              nc_inq_dimname(owner, dimensions[1], secondName) == NC_NOERR &&
              std::string(firstName) == "t" &&
              std::string(secondName) == "kl",
          "compact A0 dimensions are not [t,kl]");
  require(!findVariable(file, "Ap_real", owner, variable) &&
              !findVariable(file, "Am_real", owner, variable),
          "QG output persisted dummy Ap or Am arrays");
  for (const auto *name : {"u", "v", "eta", "pi", "psi", "qgpv",
                           "zeta_z", "ssh", "energy", "uvMax", "dye"})
    require(findVariable(file, name, owner, variable),
            std::string("QG output is missing ") + name);
  require(nc_close(file) == NC_NOERR, "unable to close QG output schema");
}

WVModelIntegratorConfiguration integratorConfiguration(
    WVModelIntegratorKind kind) {
  WVModelIntegratorConfiguration result;
  result.kind = kind;
  result.adaptive.relativeTolerance = 1e-8;
  result.adaptive.maximumStepSize = 0.005;
  result.adaptiveRK45.relativeTolerance = 1e-8;
  result.adaptiveRK45.maximumStepSize = 0.005;
  result.adaptiveRK78.relativeTolerance = 1e-8;
  result.adaptiveRK78.maximumStepSize = 0.005;
  return result;
}

void exerciseIntegrator(WVModelIntegratorKind kind, std::size_t ordinal,
                        const std::shared_ptr<const WVExtensionCatalog> &extensions,
                        const std::filesystem::path &directory) {
  const auto configuration = qgConfiguration();
  WVTransformBarotropicQGDescriptor transform;
  auto status = WVTransformBarotropicQGDescriptor::create(configuration,
                                                           transform);
  require(static_cast<bool>(status), status.message);
  const auto path = directory / ("qg-" + std::to_string(ordinal) + ".nc");
  // Persist a schedule that extends beyond the create run so the second model
  // exercises actual append output rather than state-only continuation.
  auto record = observerRecord(configuration, transform.Nkl(), path, 0.04);
  WVPortableObserverDescriptor descriptor;
  status = WVPortableObserverDescriptor::create(record, extensions, descriptor);
  require(static_cast<bool>(status), status.message);
  const auto forcing = defaultNonlinearAdvectionSchedule();
  WVModel model;
  status = WVModel::create(extensions, configuration, forcing, descriptor,
                           std::make_unique<WVReferenceFFTEngine>(),
                           integratorConfiguration(kind), model);
  require(static_cast<bool>(status), status.message);
  const auto description = descriptionFor(model.stateLayout());
  WVModelState state;
  status = WVModelState::create(
      checkpointFor(configuration, description, forcing), model.stateLayout(),
      state);
  require(static_cast<bool>(status), status.message);
  initializeAdditionalState(model, state);
  require(state.constView().coefficientFamilyCount == 1 &&
              state.constView().waveVortex.coefficients.Ap.data == nullptr &&
              state.constView().waveVortex.coefficients.Am.data == nullptr,
          "fresh model did not retain compact A0-only state");

  WVModelOutputConfiguration output;
  status = WVModelOutputConfiguration::compile(
      record, {}, {}, WVModelOutputPolicy::create, extensions, 0.0, 0.02,
      output, nullptr, false, &description);
  require(static_cast<bool>(status), status.message);
  status = model.openOutput(state, std::move(output));
  require(static_cast<bool>(status) && model.hasOutput(), status.message);
  status = model.advanceToTime(state, 0.02, 0.005);
  require(static_cast<bool>(status), status.message);
  const auto closed = model.closeOutput();
  require(static_cast<bool>(closed) && !model.hasOutput(), closed.message);
  verifyCompactSchema(path);

  WVModelOutputNetCDFInspection inspection;
  auto persistence = WVModelOutputNetCDFSink::inspect(
      {path.string()}, *extensions, inspection);
  require(static_cast<bool>(persistence), persistence.message);
  require(inspection.latestRestart.transformKind ==
                  WVPersistedTransformKind::barotropicQG &&
              inspection.latestRestart.stateDescription.coefficientFamilies
                      .size() == 1 &&
              inspection.latestRestart.stateDescription
                      .coefficientFamilies[0]
                      .identifier == "A0" &&
              inspection.latestRestart.t == 0.02,
          "QG model-output inspection lost compact restart identity");

  WVModel continued;
  WVModelState restored;
  WVModelOutputRequest request;
  request.policy = WVModelOutputPolicy::append;
  request.finalTime = 0.04;
  status = WVModel::createFromModelOutputFiles(
      extensions, {path.string()}, request,
      std::make_unique<WVReferenceFFTEngine>(), integratorConfiguration(kind),
      continued, restored);
  require(static_cast<bool>(status), status.message);
  require(relativeStateDifference(state, restored) == 0.0,
          "QG model-output restore did not reproduce the complete state");
  status = model.prepareStateAfterRestart(state);
  require(static_cast<bool>(status), status.message);
  status = continued.prepareStateAfterRestart(restored);
  require(static_cast<bool>(status), status.message);
  status = model.advanceToTime(state, 0.03, 0.005);
  require(static_cast<bool>(status), status.message);
  status = model.advanceToTime(state, 0.04, 0.005);
  require(static_cast<bool>(status), status.message);
  status = continued.advanceToTime(restored, 0.04, 0.005);
  require(static_cast<bool>(status), status.message);
  require(relativeStateDifference(state, restored) <= 1e-12,
          "segmented QG continuation differs from same-host continuation");
  const auto liveMetrics = continued.metrics(&restored);
  require(liveMetrics.outputPersistentBytes > 0 &&
              liveMetrics.outputPersistentBytes ==
                  liveMetrics.outputConfigurationPersistentBytes +
                      liveMetrics.outputEvaluationPersistentBytes +
                      liveMetrics.outputSinkPersistentBytes &&
              liveMetrics.outputEvaluation.borrowedCoefficientViewCount > 0 &&
              liveMetrics.outputEvaluation.fieldEvaluationCount > 0 &&
              liveMetrics.outputEvaluation.occurrenceWorkspaceLiveBytes == 0 &&
              liveMetrics.outputEvaluation
                      .occurrenceWorkspaceMaximumLiveBytes <=
                  liveMetrics.outputEvaluation
                      .occurrenceWorkspaceRetainedBytes,
          "QG output storage or event-scoped field reuse is not exactly "
          "accounted: total=" +
              std::to_string(liveMetrics.outputPersistentBytes) +
              " configuration=" +
              std::to_string(
                  liveMetrics.outputConfigurationPersistentBytes) +
              " evaluation=" +
              std::to_string(liveMetrics.outputEvaluationPersistentBytes) +
              " sink=" +
              std::to_string(liveMetrics.outputSinkPersistentBytes) +
              " borrowed=" +
              std::to_string(liveMetrics.outputEvaluation
                                 .borrowedCoefficientViewCount) +
              " fields=" +
              std::to_string(
                  liveMetrics.outputEvaluation.fieldEvaluationCount) +
              " reused=" +
              std::to_string(
                  liveMetrics.outputEvaluation.sharedFieldReuseCount) +
              " workspace-live=" +
              std::to_string(liveMetrics.outputEvaluation
                                 .occurrenceWorkspaceLiveBytes) +
              " workspace-max=" +
              std::to_string(liveMetrics.outputEvaluation
                                 .occurrenceWorkspaceMaximumLiveBytes) +
              " workspace-retained=" +
              std::to_string(liveMetrics.outputEvaluation
                                 .occurrenceWorkspaceRetainedBytes));
  persistence = continued.closeOutput();
  require(static_cast<bool>(persistence), persistence.message);

  const auto metrics = continued.metrics(&restored);
  require(metrics.barotropicQGKernel.planCount == 3 &&
              metrics.barotropicQGKernel.persistentFullHermitianBytes == 0 &&
              metrics.barotropicQGKernel.scratchCapacityBytes ==
                  metrics.barotropicQGKernel.halfSpectrumScratchCapacityBytes +
                      metrics.barotropicQGKernel.realScratchCapacityBytes &&
              metrics.integratedObservers.sharedRightHandSideContextCount > 0 &&
              metrics.barotropicQGForcing.physicalFieldReuseCount > 0 &&
              continued.kernelProviderIdentifier() == "reference-direct" &&
              metrics.statePersistentBytes == restored.persistentBytes(),
          "QG storage, provider, plan, or shared-context evidence is incomplete: " +
              std::to_string(metrics.barotropicQGKernel.planCount) + "," +
              std::to_string(metrics.barotropicQGKernel.persistentFullHermitianBytes) +
              "," + std::to_string(
                            metrics.integratedObservers
                                .sharedRightHandSideContextCount) +
              "," + std::to_string(
                            metrics.barotropicQGForcing
                                .physicalFieldReuseCount) +
              "," + continued.kernelProviderIdentifier());
}

void exerciseSiblingSelection(
    const std::shared_ptr<const WVExtensionCatalog> &extensions,
    const std::filesystem::path &directory) {
  const auto configuration = qgConfiguration();
  WVTransformBarotropicQGDescriptor transform;
  auto status = WVTransformBarotropicQGDescriptor::create(configuration,
                                                           transform);
  require(static_cast<bool>(status), status.message);
  const auto primaryPath = directory / "qg-sibling-primary.nc";
  const auto secondaryPath = directory / "qg-sibling-secondary.nc";
  auto record = observerRecord(configuration, transform.Nkl(), primaryPath,
                               0.03);
  record.outputFiles.front().groups = {
      {"early-restart", "early-restart", {0.01, 0.0, 0.02},
       {"coefficients", "fields", "floats", "dye"}, true},
      {"dense-fields", "dense-fields", {0.005, 0.0, 0.03}, {"fields"},
       false}};
  record.outputFiles.push_back(
      {"secondary",
       secondaryPath.string(),
       {{"latest-restart", "latest-restart", {0.005, 0.0, 0.03},
         {"coefficients", "fields", "floats", "dye"}, true},
        {"dense-fields", "dense-fields", {0.005, 0.0, 0.03}, {"fields"},
         false}}});

  WVPortableObserverDescriptor descriptor;
  status = WVPortableObserverDescriptor::create(record, extensions,
                                                 descriptor);
  require(static_cast<bool>(status), status.message);
  const auto forcing = defaultNonlinearAdvectionSchedule();
  WVModel model;
  status = WVModel::create(
      extensions, configuration, forcing, descriptor,
      std::make_unique<WVReferenceFFTEngine>(),
      integratorConfiguration(WVModelIntegratorKind::fixedRK4), model);
  require(static_cast<bool>(status), status.message);
  const auto description = descriptionFor(model.stateLayout());
  WVModelState state;
  status = WVModelState::create(
      checkpointFor(configuration, description, forcing), model.stateLayout(),
      state);
  require(static_cast<bool>(status), status.message);
  initializeAdditionalState(model, state);

  WVModelOutputConfiguration output;
  status = WVModelOutputConfiguration::compile(
      record, {}, {}, WVModelOutputPolicy::create, extensions, 0.0, 0.03,
      output, nullptr, false, &description);
  require(static_cast<bool>(status), status.message);
  status = model.openOutput(state, std::move(output));
  require(static_cast<bool>(status), status.message);
  status = model.advanceToTime(state, 0.03, 0.005);
  require(static_cast<bool>(status), status.message);
  const auto metrics = model.metrics(&state);
  require(metrics.outputDriver.maximumCoincidentRouteCount == 4 &&
              metrics.outputDriver.generatedRouteCount >
                  metrics.outputDriver.generatedEventCount &&
              metrics.outputEvaluation.occurrenceReuseCount > 0 &&
              metrics.outputEvaluation.occurrenceWorkspaceLiveBytes == 0,
          "QG sibling output did not reuse coincident event work: maximum=" +
              std::to_string(
                  metrics.outputDriver.maximumCoincidentRouteCount) +
              " routes=" +
              std::to_string(metrics.outputDriver.generatedRouteCount) +
              " events=" +
              std::to_string(metrics.outputDriver.generatedEventCount) +
              " occurrence-reuse=" +
              std::to_string(
                  metrics.outputEvaluation.occurrenceReuseCount) +
              " workspace-live=" +
              std::to_string(metrics.outputEvaluation
                                 .occurrenceWorkspaceLiveBytes));
  auto persistence = model.closeOutput();
  require(static_cast<bool>(persistence), persistence.message);

  WVModelOutputNetCDFInspection inspection;
  persistence = WVModelOutputNetCDFSink::inspect(
      {primaryPath.string(), secondaryPath.string()}, *extensions,
      inspection);
  require(static_cast<bool>(persistence), persistence.message);
  require(inspection.paths.size() == 2 &&
              inspection.observerRecord.outputFiles.size() == 2 &&
              inspection.observerRecord.outputFiles.front().groups.size() ==
                  2 &&
              inspection.observerRecord.outputFiles.back().groups.size() ==
                  2 &&
              inspection.scheduleContinuations.size() == 4 &&
              inspection.destinationProgress.size() == 4 &&
              inspection.latestRestartPath == secondaryPath.string() &&
              inspection.latestRestart.t == 0.03,
          "QG sibling inspection did not select the latest complete state");
  const auto progressFor = [&](const std::string &file,
                               const std::string &group) {
    return std::find_if(
        inspection.destinationProgress.begin(),
        inspection.destinationProgress.end(), [&](const auto &progress) {
          return progress.fileIdentifier == file &&
                 progress.groupIdentifier == group;
        });
  };
  const auto early = progressFor("primary", "early-restart");
  const auto dense = progressFor("primary", "dense-fields");
  const auto latest = progressFor("secondary", "latest-restart");
  const auto denseMirror = progressFor("secondary", "dense-fields");
  require(early != inspection.destinationProgress.end() &&
              dense != inspection.destinationProgress.end() &&
              latest != inspection.destinationProgress.end() &&
              denseMirror != inspection.destinationProgress.end() &&
              early->recordCount == 3 && dense->recordCount == 7 &&
              latest->recordCount == 7 && denseMirror->recordCount == 7 &&
              early->lastCommittedTime == 0.02 &&
              dense->lastCommittedTime == 0.03 &&
              latest->lastCommittedTime == 0.03 &&
              denseMirror->lastCommittedTime == 0.03,
          "QG sibling schedules or committed progress were not preserved");

  WVModel restoredModel;
  WVModelState restored;
  WVModelOutputRequest request;
  request.policy = WVModelOutputPolicy::append;
  request.finalTime = 0.045;
  status = WVModel::createFromModelOutputFiles(
      extensions, {primaryPath.string(), secondaryPath.string()}, request,
      std::make_unique<WVReferenceFFTEngine>(),
      integratorConfiguration(WVModelIntegratorKind::fixedRK4), restoredModel,
      restored);
  require(static_cast<bool>(status), status.message);
  require(restored.constView().coefficientFamilyCount == 1 &&
              restored.constView().waveVortex.coefficients.Ap.data == nullptr &&
              restored.constView().waveVortex.coefficients.Am.data == nullptr &&
              relativeStateDifference(state, restored) == 0.0,
          "QG sibling restore did not preserve compact complete state");
  persistence = restoredModel.closeOutput();
  require(static_cast<bool>(persistence), persistence.message);
}

void rejectUnsupportedObservers(
    const std::shared_ptr<const WVExtensionCatalog> &extensions,
    const std::filesystem::path &directory) {
  const auto configuration = qgConfiguration();
  WVTransformBarotropicQGDescriptor transform;
  auto status = WVTransformBarotropicQGDescriptor::create(configuration,
                                                           transform);
  require(static_cast<bool>(status), status.message);
  const auto requireRejected = [&](WVPortableObserverRecord record,
                                   const std::filesystem::path &path,
                                   const std::string &description) {
    WVPortableObserverDescriptor descriptor;
    auto localStatus =
        WVPortableObserverDescriptor::create(record, extensions, descriptor);
    require(static_cast<bool>(localStatus), localStatus.message);
    WVModel model;
    localStatus = WVModel::create(
        extensions, configuration, defaultNonlinearAdvectionSchedule(),
        descriptor, std::make_unique<WVReferenceFFTEngine>(),
        integratorConfiguration(WVModelIntegratorKind::fixedRK4), model);
    require(!localStatus && !std::filesystem::exists(path),
            description + " was accepted or mutated an output destination");
  };

  const auto particlePath = directory / "unsupported-qg-particles.nc";
  auto particleRecord =
      observerRecord(configuration, transform.Nkl(), particlePath, 0.02);
  auto &particles = particleRecord.observers[2];
  particles.isXYOnly = false;
  particles.z = {-10.0, -20.0};
  particles.verticalAbsoluteTolerance = 1e-5;
  particles.stateBlockIdentifiers.push_back("floats-z");
  particleRecord.stateBlocks.push_back(
      {"floats-z", WVStateScalarType::real64, {2},
       WVToleranceKind::uniformAbsolute, 1e-5,
       WVStateOwnership::integratorOwned,
       WVRestartRequirement::requiredDynamicState});
  requireRejected(std::move(particleRecord), particlePath,
                  "three-dimensional QG particles");

  const auto tracerPath = directory / "unsupported-qg-tracer.nc";
  auto tracerRecord =
      observerRecord(configuration, transform.Nkl(), tracerPath, 0.02);
  tracerRecord.observers[3].isXYOnly = false;
  tracerRecord.stateBlocks[3].dimensions = {configuration.Nx,
                                             configuration.Ny, 2};
  requireRejected(std::move(tracerRecord), tracerPath,
                  "three-dimensional QG tracer");
}

void inspectOptionalMatlabFixture(
    const std::shared_ptr<const WVExtensionCatalog> &extensions) {
  const char *path = std::getenv("WV_MATLAB_BAROTROPIC_QG_OUTPUT_FIXTURE");
  if (path == nullptr)
    return;
  WVModelOutputNetCDFInspection inspection;
  const auto status = WVModelOutputNetCDFSink::inspect({path}, *extensions,
                                                        inspection);
  require(static_cast<bool>(status), status.message);
  require(inspection.latestRestart.transformKind ==
                  WVPersistedTransformKind::barotropicQG &&
              inspection.latestRestart.stateDescription.coefficientFamilies
                      .size() == 1 &&
              inspection.latestRestart.stateDescription
                      .coefficientFamilies[0]
                      .identifier == "A0",
          "MATLAB QG fixture did not inspect as compact A0-only output");
}

} // namespace

int main() {
  try {
    TemporaryDirectory directory;
    const auto extensions = catalog();
    const std::vector<WVModelIntegratorKind> integrators = {
        WVModelIntegratorKind::fixedRK4, WVModelIntegratorKind::adaptiveRK23,
        WVModelIntegratorKind::adaptiveRK45,
        WVModelIntegratorKind::adaptiveRK78};
    for (std::size_t index = 0; index < integrators.size(); ++index)
      exerciseIntegrator(integrators[index], index, extensions,
                         directory.path);
    exerciseSiblingSelection(extensions, directory.path);
    rejectUnsupportedObservers(extensions, directory.path);
    inspectOptionalMatlabFixture(extensions);
    std::cout << "Barotropic QG model-output/restart tests passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}
