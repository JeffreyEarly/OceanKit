#pragma once

#include "WVVariableComplexBuffer.hpp"
#include <memory>
#include <limits>
#include <stdexcept>

namespace wavevortex::kernel_detail {

// Storage belongs to one kernel; validity belongs to its explicit immutable
// state evaluation. View identifiers are issued at registration, never inferred
// from coefficient addresses or time. Buffers do not alias transform scratch.
struct WVPreparedFieldKey {
    std::size_t view = 0, field = 0, component = 0;
    bool operator==(const WVPreparedFieldKey& b) const noexcept {
        return view==b.view && field==b.field && component==b.component;
    }
};

class WVPreparedFieldCache final {
public:
    struct StandaloneScope {
        WVPreparedFieldCache* cache;
        bool standalone;
        ~StandaloneScope() { if (cache && standalone) cache->clear(); }
    };
    struct Entry {
        Entry(std::size_t modalElements,std::size_t gridElements,
              WVComplexRepresentation representation)
            : data(modalElements+gridElements,representation), modalElements(modalElements), gridElements(gridElements) {}
        spectral_detail::WVVariableComplexBuffer data;
        const std::size_t modalElements,gridElements;
        WVPreparedFieldKey key{};
        bool occupied = false, modalReady = false, gridReady = false;
        WVComplexOutput modal(std::size_t offset,std::size_t count) { return data.output(offset,count); }
        WVComplexOutput grid() { return data.output(modalElements,gridElements); }
    };

    WVPreparedFieldCache(std::size_t modalElements,std::size_t gridElements,
                         WVComplexRepresentation representation)
        : modalElements_(modalElements),gridElements_(gridElements),representation_(representation) {
        if (modalElements>std::numeric_limits<std::size_t>::max()-gridElements ||
            modalElements+gridElements>std::numeric_limits<std::size_t>::max()/sizeof(WVComplex64))
            throw std::overflow_error("Prepared field extent overflow.");
        // Four prognostic fields require no allocations during a prepared RHS.
        entries_.reserve(maxEntries);
        for (std::size_t i=0;i<4;++i) entries_.push_back(makeEntry());
    }
    WVKernelStatus acquire(WVPreparedFieldKey key,Entry*& result) {
        result=nullptr;
        for (auto& entry:entries_) if (entry->occupied && entry->key==key) {
            result=entry.get(); return WVKernelStatus::ok();
        }
        Entry* available=nullptr;
        for (auto& entry:entries_) if (!entry->occupied) { available=entry.get(); break; }
        if (!available) {
            if (entries_.size()==maxEntries)
                return {WVKernelStatusCode::invalidConfiguration,"Prepared field key capacity exceeded."};
            try { entries_.push_back(makeEntry()); available=entries_.back().get(); }
            catch (const std::bad_alloc&) {
                return {WVKernelStatusCode::allocationFailure,"Prepared field allocation failed."};
            }
        }
        available->key=key; available->occupied=true;
        available->modalReady=false; available->gridReady=false;
        result=available;
        return WVKernelStatus::ok();
    }
    void invalidateView(std::size_t view) noexcept {
        for (auto& entry:entries_) if (entry->occupied && entry->key.view==view) reset(*entry);
    }
    void clear() noexcept { for (auto& entry:entries_) reset(*entry); }
    std::size_t capacityBytes() const noexcept {
        std::size_t bytes=sizeof(*this)+entries_.capacity()*sizeof(std::unique_ptr<Entry>);
        for (const auto& entry:entries_) bytes+=sizeof(Entry)+entry->data.capacityBytes();
        return bytes;
    }
private:
    // Five registered views, five requested components, nine normalized fields.
    static constexpr std::size_t maxEntries=5*5*9;
    std::unique_ptr<Entry> makeEntry() const {
        return std::make_unique<Entry>(modalElements_,gridElements_,representation_);
    }
    static void reset(Entry& entry) noexcept {
        entry.key={}; entry.occupied=false; entry.modalReady=false; entry.gridReady=false;
    }
    std::size_t modalElements_,gridElements_;
    WVComplexRepresentation representation_;
    std::vector<std::unique_ptr<Entry>> entries_;
};
} // namespace wavevortex::kernel_detail
