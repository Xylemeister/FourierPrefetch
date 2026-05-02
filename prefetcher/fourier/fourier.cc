
#include <array>
#include <unordered_map>

#include "cache.h"
#include "transform.h"

// One engine per cache instance.  A single global would be shared across
// every cache level that loads this prefetcher (L1D, L2C, LLC), mixing
// unrelated IP streams and corrupting all per-IP state.
static std::unordered_map<CACHE*, Transform::FourierPrefetchV1> engines;

void CACHE::prefetcher_initialize()
{
    engines[this];  // default-construct engine for this cache instance
}

uint32_t CACHE::prefetcher_cache_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit,
                                         bool useful_prefetch, uint8_t type,
                                         uint32_t metadata_in)
{
    auto cl_addr = addr >> LOG2_BLOCK_SIZE;
    auto& engine = engines[this];
    uint32_t meta = engine.update(cl_addr, ip, current_cycle);

    for (int64_t delta : engine.pending_deltas()) {
        if (delta == 0) continue;
        uint64_t pf_cl   = static_cast<uint64_t>(static_cast<int64_t>(cl_addr) + delta);
        uint64_t pf_addr = pf_cl << LOG2_BLOCK_SIZE;
        if ((pf_addr >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) continue;
        prefetch_line(pf_addr, (this->get_mshr_occupancy_ratio() < 0.7), meta);
    }

    return meta;
}

void CACHE::prefetcher_cycle_operate() {}

uint32_t CACHE::prefetcher_cache_fill(uint64_t addr, uint32_t set, uint32_t way,
                                      uint8_t prefetch, uint64_t evicted_addr,
                                      uint32_t metadata_in)
{
    if (prefetch)
        engines[this].record_fill(metadata_in, current_cycle);
    return metadata_in;
}

void CACHE::prefetcher_final_stats() {}
