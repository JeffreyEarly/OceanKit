#include "../../CompiledKernel/src/WVPreparedFieldCache.hpp"
#include "../../tools/compiled-kernel/tests/WVAllocationProbe.hpp"

#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>

using namespace wavevortex;
using namespace wavevortex::kernel_detail;

namespace {

void require(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

void hotEntriesArePrepared(WVComplexRepresentation representation) {
  WVPreparedFieldCache cache(12, 20, representation);
  const auto preparedBytes = cache.capacityBytes();
  std::array<WVPreparedFieldCache::Entry *, 4> entries{};

  allocationProbe::calls = 0;
  allocationProbe::counting = true;
  for (std::size_t field = 0; field < entries.size(); ++field) {
    const auto status = cache.acquire({17, field, 0}, entries[field]);
    require(bool(status) && entries[field] != nullptr,
            "A prepared prognostic entry could not be acquired.");
  }
  allocationProbe::counting = false;
  require(allocationProbe::calls == 0,
          "Acquiring the four prepared prognostic entries allocated.");
  require(cache.capacityBytes() == preparedBytes,
          "Acquiring prepared entries changed retained capacity.");

  for (std::size_t first = 0; first < entries.size(); ++first)
    for (std::size_t second = first + 1; second < entries.size(); ++second)
      require(entries[first] != entries[second],
              "Distinct prepared keys shared one cache entry.");

  entries[0]->modalReady = true;
  entries[0]->gridReady = false;
  WVPreparedFieldCache::Entry *same = nullptr;
  require(bool(cache.acquire({17, 0, 0}, same)) && same == entries[0] &&
              same->modalReady && !same->gridReady,
          "A same-key lookup lost partially published preparation state.");

  const auto modal = same->modal(0, 12);
  const auto secondModal = same->modal(6, 6);
  const auto grid = same->grid();
  require(modal.bytes == 12 * (representation == WVComplexRepresentation::interleaved
                                  ? sizeof(WVComplex64)
                                  : sizeof(double)) &&
              secondModal.bytes == 6 *
                                       (representation ==
                                                WVComplexRepresentation::interleaved
                                            ? sizeof(WVComplex64)
                                            : sizeof(double)) &&
              grid.bytes == 20 *
                                (representation == WVComplexRepresentation::interleaved
                                     ? sizeof(WVComplex64)
                                     : sizeof(double)),
          "Prepared cache views reported the wrong extent.");
  if (representation == WVComplexRepresentation::interleaved)
    require(modal.interleaved && secondModal.interleaved == modal.interleaved + 6 &&
                grid.interleaved && !modal.real && !grid.real,
            "Interleaved cache views have the wrong representation or offset.");
  else
    require(modal.real && modal.imag && secondModal.real == modal.real + 6 &&
                secondModal.imag == modal.imag + 6 && grid.real && grid.imag &&
                !modal.interleaved && !grid.interleaved,
            "Split cache views have the wrong representation or offset.");

  cache.invalidateView(99);
  require(same->occupied && same->modalReady && !same->gridReady,
          "Invalidating another view changed a prepared entry.");
  cache.invalidateView(17);
  require(!same->occupied && !same->modalReady && !same->gridReady,
          "View invalidation retained published preparation state.");

  WVPreparedFieldCache::Entry *replacement = nullptr;
  allocationProbe::calls = 0;
  allocationProbe::counting = true;
  const auto replacementStatus = cache.acquire({18, 8, 4}, replacement);
  allocationProbe::counting = false;
  require(bool(replacementStatus) && replacement != nullptr &&
              allocationProbe::calls == 0 &&
              cache.capacityBytes() == preparedBytes,
          "An invalidated prepared entry did not support allocation-free reuse.");

  replacement->modalReady = replacement->gridReady = true;
  {
    WVPreparedFieldCache::StandaloneScope nested{&cache, false};
  }
  require(replacement->occupied && replacement->modalReady && replacement->gridReady,
          "A non-standalone scope cleared an active evaluation.");
  {
    WVPreparedFieldCache::StandaloneScope standalone{&cache, true};
  }
  require(!replacement->occupied && !replacement->modalReady &&
              !replacement->gridReady,
          "A standalone scope published preparation beyond its operation.");
}

void boundedKeyTable() {
  WVPreparedFieldCache cache(1, 1, WVComplexRepresentation::interleaved);
  WVPreparedFieldCache::Entry *entry = nullptr;
  std::size_t keys = 0;
  for (std::size_t view = 0; view < 5; ++view)
    for (std::size_t component = 0; component < 5; ++component)
      for (std::size_t field = 0; field < 9; ++field) {
        require(bool(cache.acquire({view, field, component}, entry)) && entry,
                "The documented prepared-key table rejected a valid key.");
        entry->modalReady = entry->gridReady = true;
        ++keys;
      }
  require(keys == 225, "The bounded prepared-key test covered the wrong extent.");

  const auto fullBytes = cache.capacityBytes();
  allocationProbe::calls = 0;
  allocationProbe::counting = true;
  const auto overflow = cache.acquire({5, 0, 0}, entry);
  allocationProbe::counting = false;
  require(!overflow &&
              overflow.code == WVKernelStatusCode::invalidConfiguration &&
              entry == nullptr,
          "A key beyond the bounded table was accepted.");

  cache.invalidateView(2);
  allocationProbe::calls = 0;
  allocationProbe::counting = true;
  const auto reused = cache.acquire({8, 8, 4}, entry);
  allocationProbe::counting = false;
  require(bool(reused) && entry && allocationProbe::calls == 0 &&
              cache.capacityBytes() == fullBytes && !entry->modalReady &&
              !entry->gridReady,
          "Invalidated bounded storage was not reused without allocation.");
}

void allocationFailureIsRecoverable() {
  WVPreparedFieldCache cache(2, 3, WVComplexRepresentation::split);
  WVPreparedFieldCache::Entry *entry = nullptr;
  for (std::size_t field = 0; field < 4; ++field)
    require(bool(cache.acquire({1, field, 0}, entry)),
            "Prepared entry setup failed.");

  allocationProbe::failAfter = 0;
  const auto failed = cache.acquire({1, 4, 0}, entry);
  allocationProbe::failAfter = -1;
  require(!failed && failed.code == WVKernelStatusCode::allocationFailure &&
              entry == nullptr,
          "A failed extra-entry allocation was published.");
  require(bool(cache.acquire({1, 4, 0}, entry)) && entry && !entry->modalReady &&
              !entry->gridReady,
          "The cache did not recover after an allocation failure.");
}

void extentOverflowIsRejected() {
  bool threw = false;
  try {
    WVPreparedFieldCache cache(std::numeric_limits<std::size_t>::max(), 1,
                               WVComplexRepresentation::interleaved);
    (void)cache;
  } catch (const std::overflow_error &) {
    threw = true;
  }
  require(threw, "A prepared-field extent overflow was accepted.");
}

} // namespace

int main() {
  hotEntriesArePrepared(WVComplexRepresentation::interleaved);
  hotEntriesArePrepared(WVComplexRepresentation::split);
  boundedKeyTable();
  allocationFailureIsRecoverable();
  extentOverflowIsRejected();
  return 0;
}
