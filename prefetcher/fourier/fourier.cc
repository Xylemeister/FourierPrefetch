
#include <array>

#include "cache.h"
#include "transform.h"

Transform::FourierPrefetchV1 engine;

void CACHE::prefetcher_initialize() {}

uint32_t CACHE::prefetcher_cache_operate(uint64_t addr, uint64_t ip, uint8_t cache_hit,
                                         bool useful_prefetch, uint8_t type,
                                         uint32_t metadata_in)
{
    auto cl_addr = addr >> LOG2_BLOCK_SIZE;
    engine.update(cl_addr, ip);

    auto delta = engine.issue(ip);
    if (delta.has_value() && delta.value() != 0) {
        uint64_t pf_cl   = static_cast<uint64_t>(static_cast<int64_t>(cl_addr) + delta.value());
        uint64_t pf_addr = pf_cl << LOG2_BLOCK_SIZE;


        // if this works might add some info of confidence here
        prefetch_line(pf_addr, (this->get_mshr_occupancy_ratio() < 0.5), metadata_in);
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
