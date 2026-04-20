
#include <array>

#include "cache.h"
#include "transform.h"

Transform::FourierPrefetchV1 tracker;

void CACHE::prefetcher_initialize() {}

uint32_t CACHE::prefetcher_cache_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit,
                                         bool useful_prefetch, uint8_t type,
                                         uint32_t metadata_in)
{
    uint64_t cl_addr = addr >> LOG2_BLOCK_SIZE;
    bool     fill_l1 = (get_mshr_occupancy_ratio() < 0.5f);
    auto     candidates = tracker.Operate(cl_addr,ip, fill_l1);

    // Precompute all prefetch addresses (deltas are cumulative in the chain).
    uint64_t walk_cl = cl_addr;
    std::array<uint64_t, Transform::WINDOW_SIZE> pf_addrs{};
    const std::size_t n = candidates.size();
    for (std::size_t i = 0; i < n; ++i) {
        walk_cl += static_cast<uint64_t>(candidates[i].delta);
        pf_addrs[i] = walk_cl << LOG2_BLOCK_SIZE;
    }

    // Issue furthest-first so the access needing the most lead time is
    // queued first. Stop as soon as prefetch_line returns false — that
    // signals the MSHR/prefetch queue is full, so further requests would
    // be dropped anyway.
    for (std::size_t i = n; i-- > 0;)
        if (!prefetch_line(pf_addrs[i], candidates[i].fill_l1, metadata_in))
            break;

    return metadata_in;
}

void CACHE::prefetcher_cycle_operate() {}

uint32_t CACHE::prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way,
                                      uint8_t prefetch, uint64_t evicted_addr,
                                      uint32_t metadata_in)
{
    return metadata_in;
}

void CACHE::prefetcher_final_stats() {}
