#pragma once

#include "WaveVortexRuntime/WVForcing.hpp"
#include "WaveVortexRuntime/WVPortableImplementationContract.hpp"
#include "WaveVortexKernel/WVTransformBarotropicQGKernel.hpp"
#include "WaveVortexKernel/WVStratifiedModalSource.hpp"

#include "WaveVortexKernel/WVTransformHydrostaticKernel.hpp"
#include "WaveVortexKernel/WVTransformBoussinesqKernel.hpp"

#include <cstdint>
#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace wavevortex::runtime {

class WVBarotropicQGForcing;
class WVStratifiedQGForcing;

enum class WVForcingPersistenceEncoding : std::uint8_t {
  realVariable,
  logicalVariable,
  textAttribute,
  zeroBasedIndexVariable,
  complexVariable
};

enum class WVForcingDimensionRule : std::uint8_t {
  scalar,
  ownLength,
  referencedLength,
  horizontalYX,
  componentPair
};

struct WVForcingPersistenceField {
  WVForcingPersistenceEncoding encoding =
      WVForcingPersistenceEncoding::realVariable;
  std::string recordName;
  std::string imaginaryRecordName;
  std::string netcdfName;
  WVForcingDimensionRule dimensions = WVForcingDimensionRule::scalar;
  std::string dimensionReference;
  bool optional = false;
  bool nonnegative = false;
  bool positive = false;
  bool allowInfinity = false;
};

struct WVForcingPersistenceSchema {
  std::vector<WVForcingPersistenceField> fields;
  bool writesNameAttribute = false;
};

using WVForcingFactory = std::function<WVKernelStatus(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &,
    const WVForcingPreparation &, std::unique_ptr<WVForcing> &)>;

using WVBarotropicQGForcingPreflight = std::function<WVKernelStatus(
    const WVFrozenForcingEntry &, std::size_t)>;

using WVBarotropicQGForcingFactory = std::function<WVKernelStatus(
    const WVFrozenForcingEntry &,
    const WVTransformBarotropicQGDescriptor &, bool,
    std::unique_ptr<WVBarotropicQGForcing> &)>;

using WVStratifiedQGForcingFactory = std::function<WVKernelStatus(
    const WVFrozenForcingEntry &, const WVStratifiedModalGeometry &, bool,
    std::unique_ptr<WVStratifiedQGForcing> &)>;

using WVHydrostaticForcingFactory = std::function<WVKernelStatus(const WVFrozenForcingEntry&, WVTransformHydrostaticKernel&, const WVForcingPreparation&, std::unique_ptr<WVForcing>&)>;
using WVBoussinesqForcingFactory = std::function<WVKernelStatus(const WVFrozenForcingEntry&, WVTransformBoussinesqKernel&, const WVForcingPreparation&, std::unique_ptr<WVForcing>&)>;

// Data-only evaluation requirements declared by each forcing registration.
// Constant Laplacian entries are [horizontal, vertical]. Variable-
// stratification grid-calculus entries are [field*4+kind], with fields
// [u,v,w,eta] and kinds [Dxx,Dyy,F-Dzz,G-Dzz].
struct WVForcingEvaluationDependencies {
  std::size_t nonlinearUseCount = 0;
  bool horizontalMaximumNeeded = false;
  std::array<std::size_t, 2> constantLaplacianUseCount{};
  std::array<std::size_t, 16> hydrostaticGridCalculusUseCount{};
  std::array<std::size_t, 16> boussinesqGridCalculusUseCount{};
};

struct WVForcingFactoryRegistration {
  std::string matlabClassName;
  std::uint32_t contractVersion = WVPortablePairContractVersion;
  std::vector<std::string> forcingTypes;
  std::string defaultName;
  WVForcingStage stage = WVForcingStage::spatial;
  std::uint8_t priority = 255;
  WVForcingPersistenceSchema persistence;
  WVForcingFactory factory;
  bool isSupported = true;
  std::string unavailabilityReason;
  bool providesAdaptiveDamping = false;
  WVForcingStage barotropicQGStage = WVForcingStage::spatial;
  WVBarotropicQGForcingPreflight barotropicQGPreflight;
  WVBarotropicQGForcingFactory barotropicQGFactory;
  // Called during model preflight, before FFT plans or mutable output.
  std::function<WVKernelStatus(const WVFrozenForcingEntry &, bool)> modelPreflight = {};
  std::function<WVKernelStatus(const WVFrozenForcingEntry &,
      const WVTransformConstantStratificationDescriptor &, WVForcingPreparation &)>
      prepareConstantResolution = {};
  WVForcingStage stratifiedQGStage = WVForcingStage::spatial;
  WVBarotropicQGForcingPreflight stratifiedQGPreflight = {};
  WVStratifiedQGForcingFactory stratifiedQGFactory = {};
  WVHydrostaticForcingFactory hydrostaticFactory = {};
  std::function<WVKernelStatus(const WVFrozenForcingEntry &,
      const WVTransformHydrostaticKernel &, WVForcingPreparation &)>
      prepareHydrostaticResolution = {};
  WVBoussinesqForcingFactory boussinesqFactory = {};
  std::function<WVKernelStatus(const WVFrozenForcingEntry &,
      const WVTransformBoussinesqKernel &, WVForcingPreparation &)>
      prepareBoussinesqResolution = {};
  // Data-only diagnostic metadata may be resolved before creating FFT plans.
  bool supportsTendencyDiagnostics = false;
  WVForcingEvaluationDependencies evaluationDependencies = {};
};

std::vector<WVForcingFactoryRegistration> builtInForcingFactories();

WVFrozenForcingSchedule defaultNonlinearAdvectionSchedule();

} // namespace wavevortex::runtime
