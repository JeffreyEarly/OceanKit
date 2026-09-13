#pragma once

#include "WaveVortexRuntime/WVForcingContracts.hpp"

namespace wavevortex::runtime::detail {

WVKernelStatus preflightLaplacianDamping(const WVFrozenForcingEntry &, bool);
WVKernelStatus preflightVerticalDiffusivity(const WVFrozenForcingEntry &, bool);
WVKernelStatus createHorizontalDamping(const WVFrozenForcingEntry &, const WVTransformConstantStratificationDescriptor &, const WVForcingPreparation &, std::unique_ptr<WVForcing> &);
WVKernelStatus createVerticalDamping(const WVFrozenForcingEntry &, const WVTransformConstantStratificationDescriptor &, const WVForcingPreparation &, std::unique_ptr<WVForcing> &);
WVKernelStatus createVerticalDiffusivity(const WVFrozenForcingEntry &, const WVTransformConstantStratificationDescriptor &, const WVForcingPreparation &, std::unique_ptr<WVForcing> &);
WVKernelStatus prepareExplicitAntialiasing(const WVFrozenForcingEntry &, const WVTransformConstantStratificationDescriptor &, WVForcingPreparation &);
WVKernelStatus preflightExplicitAntialiasing(const WVFrozenForcingEntry &, bool);
WVKernelStatus createExplicitAntialiasing(
    const WVFrozenForcingEntry &, const WVTransformConstantStratificationDescriptor &,
    const WVForcingPreparation &, std::unique_ptr<WVForcing> &);
WVKernelStatus createBarotropicQGExplicitAntialiasing(
    const WVFrozenForcingEntry &, const WVTransformBarotropicQGDescriptor &,
    bool, std::unique_ptr<WVBarotropicQGForcing> &);
WVKernelStatus preflightBarotropicQGExplicitAntialiasing(const WVFrozenForcingEntry &, std::size_t);

WVKernelStatus createNonlinearAdvectionForcing(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &, const WVForcingPreparation &,
    std::unique_ptr<WVForcing> &);
WVKernelStatus createAdaptiveDampingForcing(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &, const WVForcingPreparation &,
    std::unique_ptr<WVForcing> &);
WVKernelStatus createFixedAmplitudeForcing(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &, const WVForcingPreparation &,
    std::unique_ptr<WVForcing> &);
WVKernelStatus createQuadraticBottomFriction(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &, const WVForcingPreparation &,
    std::unique_ptr<WVForcing> &);
WVKernelStatus createLinearBottomFriction(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &, const WVForcingPreparation &,
    std::unique_ptr<WVForcing> &);
WVKernelStatus createPseudoTopographicForcing(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &, const WVForcingPreparation &,
    std::unique_ptr<WVForcing> &);
WVKernelStatus createBetaPlaneForcing(
    const WVFrozenForcingEntry &,
    const WVTransformConstantStratificationDescriptor &, const WVForcingPreparation &,
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
WVKernelStatus preflightStratifiedQGEmptyForcing(
    const WVFrozenForcingEntry &, std::size_t);
WVKernelStatus preflightStratifiedQGFixedAmplitude(
    const WVFrozenForcingEntry &, std::size_t);
WVKernelStatus preflightStratifiedQGScalarForcing(
    const WVFrozenForcingEntry &, std::size_t);
WVKernelStatus createStratifiedQGNonlinearAdvection(
    const WVFrozenForcingEntry &,
    const WVStratifiedModalGeometry &, bool,
    std::unique_ptr<WVStratifiedQGForcing> &);
WVKernelStatus createStratifiedQGAdaptiveDamping(
    const WVFrozenForcingEntry &,
    const WVStratifiedModalGeometry &, bool,
    std::unique_ptr<WVStratifiedQGForcing> &);
WVKernelStatus createStratifiedQGFixedAmplitude(
    const WVFrozenForcingEntry &,
    const WVStratifiedModalGeometry &, bool,
    std::unique_ptr<WVStratifiedQGForcing> &);
WVKernelStatus createStratifiedQGQuadraticBottomFriction(
    const WVFrozenForcingEntry &,
    const WVStratifiedModalGeometry &, bool,
    std::unique_ptr<WVStratifiedQGForcing> &);
WVKernelStatus createStratifiedQGLinearBottomFriction(
    const WVFrozenForcingEntry &,
    const WVStratifiedModalGeometry &, bool,
    std::unique_ptr<WVStratifiedQGForcing> &);
WVKernelStatus createStratifiedQGBetaPlanePVAdvection(
    const WVFrozenForcingEntry &,
    const WVStratifiedModalGeometry &, bool,
    std::unique_ptr<WVStratifiedQGForcing> &);

WVKernelStatus preflightStratifiedQGExplicitAntialiasing(const WVFrozenForcingEntry &,std::size_t);
WVKernelStatus preflightStratifiedQGVerticalDiffusivity(const WVFrozenForcingEntry &,std::size_t);
WVKernelStatus createStratifiedQGExplicitAntialiasing(const WVFrozenForcingEntry &,const WVStratifiedModalGeometry &,bool,std::unique_ptr<WVStratifiedQGForcing> &);
WVKernelStatus createStratifiedQGVerticalDiffusivity(const WVFrozenForcingEntry &,const WVStratifiedModalGeometry &,bool,std::unique_ptr<WVStratifiedQGForcing> &);

WVKernelStatus createHydrostaticHorizontalDamping(const WVFrozenForcingEntry&,WVTransformHydrostaticKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createHydrostaticVerticalDamping(const WVFrozenForcingEntry&,WVTransformHydrostaticKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createHydrostaticVerticalDiffusivity(const WVFrozenForcingEntry&,WVTransformHydrostaticKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createHydrostaticNonlinearAdvectionForcing(const WVFrozenForcingEntry&,WVTransformHydrostaticKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createHydrostaticAdaptiveDampingForcing(const WVFrozenForcingEntry&,WVTransformHydrostaticKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createHydrostaticQuadraticBottomFriction(const WVFrozenForcingEntry&,WVTransformHydrostaticKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createHydrostaticLinearBottomFriction(const WVFrozenForcingEntry&,WVTransformHydrostaticKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createHydrostaticFixedAmplitudeForcing(const WVFrozenForcingEntry&,WVTransformHydrostaticKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createHydrostaticBetaPlaneForcing(const WVFrozenForcingEntry&,WVTransformHydrostaticKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createHydrostaticPseudoTopographicForcing(const WVFrozenForcingEntry&,WVTransformHydrostaticKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus prepareHydrostaticExplicitAntialiasing(const WVFrozenForcingEntry&, const WVTransformHydrostaticKernel&, WVForcingPreparation&);
WVKernelStatus createHydrostaticExplicitAntialiasing(const WVFrozenForcingEntry&,WVTransformHydrostaticKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createBoussinesqHorizontalDamping(const WVFrozenForcingEntry&,WVTransformBoussinesqKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createBoussinesqVerticalDamping(const WVFrozenForcingEntry&,WVTransformBoussinesqKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createBoussinesqVerticalDiffusivity(const WVFrozenForcingEntry&,WVTransformBoussinesqKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createBoussinesqNonlinearAdvectionForcing(const WVFrozenForcingEntry&,WVTransformBoussinesqKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createBoussinesqAdaptiveDampingForcing(const WVFrozenForcingEntry&,WVTransformBoussinesqKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createBoussinesqQuadraticBottomFriction(const WVFrozenForcingEntry&,WVTransformBoussinesqKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createBoussinesqLinearBottomFriction(const WVFrozenForcingEntry&,WVTransformBoussinesqKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createBoussinesqFixedAmplitudeForcing(const WVFrozenForcingEntry&,WVTransformBoussinesqKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createBoussinesqBetaPlaneForcing(const WVFrozenForcingEntry&,WVTransformBoussinesqKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus createBoussinesqPseudoTopographicForcing(const WVFrozenForcingEntry&,WVTransformBoussinesqKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);
WVKernelStatus prepareBoussinesqExplicitAntialiasing(const WVFrozenForcingEntry&, const WVTransformBoussinesqKernel&, WVForcingPreparation&);
WVKernelStatus createBoussinesqExplicitAntialiasing(const WVFrozenForcingEntry&,WVTransformBoussinesqKernel&,const WVForcingPreparation&,std::unique_ptr<WVForcing>&);

} // namespace wavevortex::runtime::detail
