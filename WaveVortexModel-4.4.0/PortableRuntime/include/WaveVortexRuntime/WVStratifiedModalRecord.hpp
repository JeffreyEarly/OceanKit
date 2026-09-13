#pragma once

#include "WaveVortexRuntime/WVCheckpointReader.hpp"
#include "WaveVortexKernel/WVOwnedStratifiedModalSource.hpp"
#include <memory>

namespace wavevortex::runtime {

// Scientific record semantics, separate from checkpoint and execution contracts.
inline constexpr const char* WVStratifiedModalRecordContract = "wave-vortex-stratified-modal-record-v1";

class WVStratifiedModalRecord final : public WVStratifiedModalSource {
public:
    // Scientific ownership is fixed once published to consumers.
    WVStratifiedModalRecord(const WVStratifiedModalRecord&) = delete;
    WVStratifiedModalRecord& operator=(const WVStratifiedModalRecord&) = delete;
    WVStratifiedModalRecord(WVStratifiedModalRecord&&) = delete;
    WVStratifiedModalRecord& operator=(WVStratifiedModalRecord&&) = delete;
    const WVStratifiedModalGeometry& geometry() const noexcept { return source_->geometry(); }
    // Opaque MATLAB .mat payload, retained only for round-trip persistence.
    // Portable numerical execution never deserializes or executes this code.
    const std::vector<unsigned char>& N2FunctionPayload() const noexcept { return N2FunctionPayload_; }
    const std::vector<double>& PF0inv() const noexcept { return source_->PF0inv(); }
    const std::vector<double>& QG0inv() const noexcept { return source_->QG0inv(); }
    const std::vector<double>& PF0() const noexcept { return source_->PF0(); }
    const std::vector<double>& QG0() const noexcept { return source_->QG0(); }
    const std::vector<double>& PFpmInv() const noexcept { return source_->PFpmInv(); }
    const std::vector<double>& QGpmInv() const noexcept { return source_->QGpmInv(); }
    const std::vector<double>& PFpm() const noexcept { return source_->PFpm(); }
    const std::vector<double>& QGpm() const noexcept { return source_->QGpm(); }
    const std::vector<double>& QGwg() const noexcept { return source_->QGwg(); }
    // Exact identity of this immutable in-process scientific record. No reuse by
    // path, radius or approximate matrix values; rereading creates a new record.
    const std::string& sourceIdentity() const noexcept { return source_->sourceIdentity(); }
    const std::string& modeSetIdentity() const noexcept { return source_->modeSetIdentity(); }
    std::size_t persistentBytes() const noexcept;
    // Exact persisted wave-group membership for Boussinesq; one shared group
    // for Hydrostatic/SQG. Balanced operators always cover all columns.
    const std::vector<WVScientificModalGroup>& groups() const noexcept { return source_->groups(); }
    // Borrowed const view of persisted preconditioned values, valid while this
    // scientific record lives. No normalization conversion or derived product.
    WVKernelStatus matrixView(WVStratifiedScientificMatrix matrix, WVVerticalMatrixRecord& result) const;

    WVKernelStatus prepareVertical(WVStratifiedModalOperator operation,
        WVComplexLayout input, WVComplexLayout output,
        std::unique_ptr<WVVerticalMatrixBackend> backend,
        std::unique_ptr<WVPreparedVerticalOperator>& result,
        WVAccumulation accumulation = WVAccumulation::overwrite) const;
    WVKernelStatus horizontalSpecification(std::size_t planes,
        WVComplexRepresentation representation, const std::string& family,
        WVRetainedHorizontalSpecification& result) const;
private:
    friend class WVStratifiedModalReader;
    WVStratifiedModalRecord() = default;
    std::vector<unsigned char> N2FunctionPayload_;
    std::shared_ptr<const WVOwnedStratifiedModalSource> source_;
};

class WVStratifiedModalReader {
public:
    // Reads geometry/modal data only. Does not promise transform execution or
    // decode forcings/coefficient state; those remain checkpoint/kernel concerns.
    // Both calls validate all matrix payloads. Failure preserves existing output.
    static WVCheckpointStatus inspect(const std::string& path, WVStratifiedModalInspection& result);
    static WVCheckpointStatus read(const std::string& path, std::shared_ptr<const WVStratifiedModalRecord>& result);
};

} // namespace wavevortex::runtime
