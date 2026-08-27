#pragma once

#include "WaveVortexRuntime/WVForcingContracts.hpp"

namespace wavevortex::runtime::detail {

WVKernelStatus createNonlinearAdvectionForcing(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &, bool,
    std::unique_ptr<WVForcing> &);
WVKernelStatus createAdaptiveDampingForcing(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &, bool,
    std::unique_ptr<WVForcing> &);
WVKernelStatus createFixedAmplitudeForcing(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &, bool,
    std::unique_ptr<WVForcing> &);
WVKernelStatus createQuadraticBottomFriction(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &, bool,
    std::unique_ptr<WVForcing> &);
WVKernelStatus createLinearBottomFriction(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &, bool,
    std::unique_ptr<WVForcing> &);
WVKernelStatus createPseudoTopographicForcing(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &, bool,
    std::unique_ptr<WVForcing> &);
WVKernelStatus createBetaPlaneForcing(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &, bool,
    std::unique_ptr<WVForcing> &);

WVKernelStatus preflightBarotropicQGEmptyForcing(
    const WVFrozenForcingEntry &, std::size_t);
WVKernelStatus preflightBarotropicQGFixedAmplitude(
    const WVFrozenForcingEntry &, std::size_t);
WVKernelStatus preflightBarotropicQGScalarForcing(
    const WVFrozenForcingEntry &, std::size_t);
WVKernelStatus createBarotropicQGNonlinearAdvection(
    const WVFrozenForcingEntry &,
    const WVTransformBarotropicQGDescriptor &, bool,
    std::unique_ptr<WVBarotropicQGForcing> &);
WVKernelStatus createBarotropicQGAdaptiveDamping(
    const WVFrozenForcingEntry &,
    const WVTransformBarotropicQGDescriptor &, bool,
    std::unique_ptr<WVBarotropicQGForcing> &);
WVKernelStatus createBarotropicQGFixedAmplitude(
    const WVFrozenForcingEntry &,
    const WVTransformBarotropicQGDescriptor &, bool,
    std::unique_ptr<WVBarotropicQGForcing> &);
WVKernelStatus createBarotropicQGQuadraticBottomFriction(
    const WVFrozenForcingEntry &,
    const WVTransformBarotropicQGDescriptor &, bool,
    std::unique_ptr<WVBarotropicQGForcing> &);
WVKernelStatus createBarotropicQGLinearBottomFriction(
    const WVFrozenForcingEntry &,
    const WVTransformBarotropicQGDescriptor &, bool,
    std::unique_ptr<WVBarotropicQGForcing> &);
WVKernelStatus createBarotropicQGBetaPlanePVAdvection(
    const WVFrozenForcingEntry &,
    const WVTransformBarotropicQGDescriptor &, bool,
    std::unique_ptr<WVBarotropicQGForcing> &);

} // namespace wavevortex::runtime::detail
