#pragma once
#include <array>
#include "WVForcingTendency.hpp"
#include "WaveVortexRuntime/WVVariableEvaluation.hpp"
#include "WVForcingEngine.hpp"
#include "WaveVortexKernel/WVTransformHydrostaticKernel.hpp"
#include "WaveVortexRuntime/WVVariableKernelServices.hpp"
namespace wavevortex::runtime {
namespace detail { class WVForcingDiagnosticBinding; }
// Reuses the resolved WVForcing services; scientific operators are supplied by
// the immutable Hydrostatic source and prepared once at construction.
class WVHydrostaticForcingEngine final {
public:
    // Indices refer to the already resolved stage/priority order.
    std::size_t forcingCount() const noexcept { return forcing_.size(); }
    const WVForcing* forcingInstance(std::size_t index) const noexcept {
        return index<forcing_.size() ? forcing_[index].get() : nullptr;
    }
    const WVForcingEvaluationDependencies*
    forcingEvaluationDependencies(std::size_t index) const noexcept;
    // Optional u/v/w/eta fields must describe this exact state and time.
    // They are borrowed for this invocation and must not alias state or outputs.
    WVKernelStatus evaluateForcingTendencies(const WVState&,
        const WVForcingTendencyOutput*,std::size_t,
        const WVRealFieldBundleConstView* preparedPhysical = nullptr,
        detail::WVForcingDiagnosticWorkspace* session = nullptr);
    const WVForcingTendencyMetrics& tendencyMetrics() const noexcept { return tendencyMetrics_; }

    static WVKernelStatus validateSchedule(const WVStratifiedModalGeometry&,const WVFrozenForcingSchedule&,WVShape2D,const WVExtensionCatalog&);
    static WVKernelStatus create(std::shared_ptr<const WVStratifiedModalSource>,const WVFrozenForcingSchedule&,std::shared_ptr<const WVExtensionCatalog>,std::unique_ptr<WVFFTEngine>,std::unique_ptr<WVHydrostaticForcingEngine>&,const WVVariableKernelServices& services = {});
    // Borrowed coefficients must stay immutable until endStateEvaluation().
    WVKernelStatus beginStateEvaluation(const WVState&);
    void endStateEvaluation() noexcept;
    bool stateEvaluationActive() const noexcept { return evaluation_.active(); }
    WVKernelStatus validateStateEvaluation(const WVState&) const;
    WVKernelStatus setVariableEvaluationPolicy(WVVariableEvaluationPolicy policy);
    WVKernelStatus validateVariableEvaluationPolicyChange(
        WVVariableEvaluationPolicy policy) const noexcept;
    const WVVariableEvaluationMetrics& variableEvaluationMetrics() const noexcept { return evaluation_.metrics(); }
    WVKernelStatus nonlinearFlux(const WVState&,WVFlux&);
    WVKernelStatus physicalFields(const WVState&,WVRealFieldBundleConstView&);
    WVStateConstraintResult restoreForcingAmplitudes(WVMutableCoefficients&);
    WVKernelStatus createErrorPolicy(double,std::unique_ptr<WVIntegrationErrorPolicy>&) const;
    WVTransformHydrostaticKernel& kernel() noexcept { return *kernel_; }
    const WVTransformHydrostaticKernel& kernel() const noexcept { return *kernel_; }
    const WVForcingPreparation& preparation() const noexcept { return preparation_; }
    const WVForcingEngineMetrics& metrics() const noexcept { return metrics_; }
    const std::string& scheduleIdentifier() const noexcept { return scheduleIdentifier_; }
    std::size_t persistentBytes() const noexcept;
    WVKernelStatus speedMaxima(const WVState&,double& uv,double& w);
    // Linear evolution retains instances for diagnostics and amplitude constraints.
    // Only their ordinary coefficient RHS contributions are disabled.
    void setLinearDynamics(bool linear) noexcept { linearDynamics_ = linear; }

private:
    friend class detail::WVForcingDiagnosticBinding;
    WVKernelStatus evaluateForcingTendenciesImpl(const WVState&,
        const WVForcingTendencyOutput*,std::size_t,
        const WVRealFieldBundleConstView*,detail::WVForcingDiagnosticWorkspace*,
        WVFlux*);
    bool linearDynamics_ = false;
    WVHydrostaticForcingEngine()=default;
    WVKernelStatus initialize(const WVFrozenForcingSchedule&);
    WVKernelStatus addNonlinearFlux(const WVState&,WVFlux&);
    WVKernelStatus gridSecondDerivative(WVRealVolumeConstView,std::size_t field,
        std::size_t kind,WVHydrostaticFamily,const double*& result);
    WVKernelStatus addProjectedSpatialTendency(const WVState&,WVRealFieldBundleConstView,WVFlux&);
    WVRealFieldBundleView clearedSpatialTendency();
    WVKernelStatus addLaplacianDamping(const WVState&,double,double,WVLaplacianDirection,WVFlux&);
    WVKernelStatus addVerticalDiffusivity(const WVState&,double,bool,WVFlux&);
    WVKernelStatus addAdaptiveDamping(const WVState&,const std::vector<double>&,WVFlux&);
    WVKernelStatus addPseudoTopographicGeneration(const WVState&,const WVPseudoTopographicOperators&,WVFlux&);
    std::unique_ptr<WVTransformHydrostaticKernel> kernel_;
    std::shared_ptr<const WVExtensionCatalog> catalog_;
    std::vector<std::unique_ptr<WVForcing>> forcing_;
    std::vector<double> physical_,spatial_,derivative_;
    std::vector<double> gridCalculus_;
    std::array<std::size_t,16> gridCalculusUseCount_{};
    std::array<int,16> gridCalculusSlots_{};
    std::vector<WVComplex64> temporary_, nonlinearCache_;
    WVForcingPreparation preparation_;
    WVForcingEngineMetrics metrics_;
    std::string scheduleIdentifier_;
    bool executing_=false;
    WVVariableEvaluationPolicy evaluationPolicy_=WVVariableEvaluationPolicy::reuse;
    WVVariableEvaluationContext evaluation_;
    std::vector<std::pair<WVVariableEvaluationKey,std::size_t>> physicalGroup_;
    WVState evaluationState_{};
    double horizontalMaximum_=0,verticalMaximum_=0;
    WVKernelStatus ensurePhysicalField(const WVState&,std::size_t);
    WVKernelStatus horizontalSpeedMaximum(const WVState&,double&);
    detail::WVForcingDiagnosticWorkspace* diagnosticWorkspace_=nullptr;
    WVForcingTendencyMetrics tendencyMetrics_;
    friend class WVForcingExecutionContext;
};
}
