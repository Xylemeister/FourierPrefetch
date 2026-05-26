#include <cstdint>
#include <unordered_map>

#include "cache.h"
#include "msl/lru_table.h"
#include "transform.h"

namespace {

constexpr std::size_t TRACKER_SETS = 256;
constexpr std::size_t TRACKER_WAYS = 4;

constexpr double MSHR_DEMOTE_RATIO = 0.5;

class PeriodicityPrefetcher
{
    struct Entry {
        uint64_t                                              ip            = 0;
        uint64_t                                              last_cl_addr  = 0;
        periodicity::PeriodicityBuf<periodicity::WINDOW_SIZE> buf{};
        uint16_t                                              frontier_step = 0;

        auto index() const { return ip; }
        auto tag()   const { return ip; }
    };

    champsim::msl::lru_table<Entry> table_{TRACKER_SETS, TRACKER_WAYS};

public:
    template <typename EmitFn>
    void Update(uint64_t cl_addr, uint64_t ip, uint64_t /*demand_page*/,
                EmitFn&& /*emit*/)
    {
        Entry probe;
        probe.ip = ip;
        auto found = table_.check_hit(probe);

        if (!found) {
            Entry seed;
            seed.ip           = ip;
            seed.last_cl_addr = cl_addr;
            table_.fill(seed);
            return;
        }

        AdvanceStream(*found, cl_addr);
        if (found->frontier_step > 0) --found->frontier_step;

        table_.fill(*found);
    }

private:
    static void AdvanceStream(Entry& e, uint64_t cl_addr)
    {
        const int64_t stride = static_cast<int64_t>(cl_addr)
                             - static_cast<int64_t>(e.last_cl_addr);
        e.buf.Insert(stride);
        e.last_cl_addr = cl_addr;
    }
};

std::unordered_map<CACHE*, PeriodicityPrefetcher> engines;

}  // namespace

void CACHE::prefetcher_initialize()
{
    engines[this];
}

uint32_t CACHE::prefetcher_cache_operate(uint64_t addr, uint64_t ip,
                                         uint8_t /*cache_hit*/,
                                         bool /*useful_prefetch*/,
                                         uint8_t /*type*/,
                                         uint32_t metadata_in)
{
    const uint64_t cl_addr     = addr >> LOG2_BLOCK_SIZE;
    const uint64_t demand_page = addr >> LOG2_PAGE_SIZE;

    auto& engine = engines[this];
    engine.Update(cl_addr, ip, demand_page,
        [this](uint64_t pf_addr, bool l2_only) -> bool {
            const bool fill_here = !l2_only
                && this->get_mshr_occupancy_ratio() < MSHR_DEMOTE_RATIO;
            return prefetch_line(pf_addr, fill_here, 0);
        });

    return metadata_in;
}

void CACHE::prefetcher_cycle_operate() {}

uint32_t CACHE::prefetcher_cache_fill(uint64_t /*addr*/, uint32_t /*set*/,
                                      uint32_t /*way*/, uint8_t /*prefetch*/,
                                      uint64_t /*evicted_addr*/,
                                      uint32_t metadata_in)
{
    return metadata_in;
}

void CACHE::prefetcher_final_stats() {}
