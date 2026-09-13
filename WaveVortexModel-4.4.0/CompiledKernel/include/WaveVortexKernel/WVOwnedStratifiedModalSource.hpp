#pragma once

#include "WVStratifiedModalSource.hpp"

#include <memory>

namespace wavevortex {

// Complete, already-solved stratified modal state. Matrix storage is MATLAB
// compatible column-major storage. Shared matrices are [Nz,Nj] inverses and
// [Nj,Nz] projections. Wave matrices append a K2unique group dimension.
struct WVStratifiedModalArrays {
    WVStratifiedModalGeometry geometry;
    std::vector<double> PF0inv, QG0inv, PF0, QG0;
    std::vector<double> PFpmInv, QGpmInv, PFpm, QGpm, QGwg;
};

// NetCDF-independent immutable owner for modal arrays supplied by an
// application that has already solved the vertical eigenproblems.
class WVOwnedStratifiedModalSource final : public WVStratifiedModalSource {
public:
    WVOwnedStratifiedModalSource(const WVOwnedStratifiedModalSource&) = delete;
    WVOwnedStratifiedModalSource& operator=(const WVOwnedStratifiedModalSource&) = delete;
    WVOwnedStratifiedModalSource(WVOwnedStratifiedModalSource&&) = delete;
    WVOwnedStratifiedModalSource& operator=(WVOwnedStratifiedModalSource&&) = delete;

    static WVKernelStatus create(WVStratifiedModalArrays arrays,
        std::shared_ptr<const WVOwnedStratifiedModalSource>& result);

    const WVStratifiedModalGeometry& geometry() const noexcept override { return arrays_.geometry; }
    const std::string& sourceIdentity() const noexcept override { return sourceIdentity_; }
    const std::string& modeSetIdentity() const noexcept override { return modeSetIdentity_; }
    std::size_t persistentBytes() const noexcept override;

    const std::vector<double>& PF0inv() const noexcept { return arrays_.PF0inv; }
    const std::vector<double>& QG0inv() const noexcept { return arrays_.QG0inv; }
    const std::vector<double>& PF0() const noexcept { return arrays_.PF0; }
    const std::vector<double>& QG0() const noexcept { return arrays_.QG0; }
    const std::vector<double>& PFpmInv() const noexcept { return arrays_.PFpmInv; }
    const std::vector<double>& QGpmInv() const noexcept { return arrays_.QGpmInv; }
    const std::vector<double>& PFpm() const noexcept { return arrays_.PFpm; }
    const std::vector<double>& QGpm() const noexcept { return arrays_.QGpm; }
    const std::vector<double>& QGwg() const noexcept { return arrays_.QGwg; }
    const std::vector<WVScientificModalGroup>& groups() const noexcept { return groups_; }

    WVKernelStatus matrixView(WVStratifiedScientificMatrix matrix,
        WVVerticalMatrixRecord& result) const;
    WVKernelStatus prepareVertical(WVStratifiedModalOperator operation,
        WVComplexLayout input, WVComplexLayout output,
        std::unique_ptr<WVVerticalMatrixBackend> backend,
        std::unique_ptr<WVPreparedVerticalOperator>& result,
        WVAccumulation accumulation = WVAccumulation::overwrite) const override;
    WVKernelStatus horizontalSpecification(std::size_t planes,
        WVComplexRepresentation representation, const std::string& family,
        WVRetainedHorizontalSpecification& result) const override;

private:
    WVOwnedStratifiedModalSource() = default;
    WVKernelStatus prepareWaveVertical(WVStratifiedModalOperator operation,
        WVComplexLayout input, WVComplexLayout output,
        std::unique_ptr<WVVerticalMatrixBackend> backend,
        std::unique_ptr<WVPreparedVerticalOperator>& result,
        WVAccumulation accumulation) const;

    WVStratifiedModalArrays arrays_;
    std::string sourceIdentity_, modeSetIdentity_;
    std::vector<WVScientificModalGroup> groups_;
};

} // namespace wavevortex
