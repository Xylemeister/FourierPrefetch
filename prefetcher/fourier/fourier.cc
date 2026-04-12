
#include <cassert>

#include "cache.h"
#include "transform.h"

Transform::FourierPrefetchV1 tracker;

void CACHE::prefetcher_initialize() {}

uint32_t CACHE::prefetcher_cache_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit,
                                         bool useful_prefetch, uint8_t type,
                                         uint32_t metadata_in)
{
    uint64_t cl_addr   = addr >> LOG2_BLOCK_SIZE;
    auto     candidates = tracker.Operate(cl_addr, ip);

    // Walk the candidates, accumulating a running address offset.
    // Each candidate's delta is relative to the previous address in the chain.
    uint64_t walk_cl = cl_addr;
    for (auto& c : candidates) {
        walk_cl += static_cast<uint64_t>(c.delta);

        uint64_t pf_addr = walk_cl << LOG2_BLOCK_SIZE;

        // Do not cross virtual page boundaries
        if ((pf_addr >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) continue;

        prefetch_line(pf_addr, c.fill_l1, metadata_in);
    }

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
