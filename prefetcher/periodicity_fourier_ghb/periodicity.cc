#include <array>
#include <cstdint>
#include <unordered_map>

#include "cache.h"
#include "transform.h"

namespace {

constexpr std::size_t GHB_SIZE = 512;
constexpr std::size_t IT_SIZE  = 512;
constexpr std::size_t IT_MASK  = IT_SIZE - 1;

constexpr double MSHR_DEMOTE_RATIO = 0.5;

class PeriodicityPrefetcher
{
    struct GhbEntry {
        uint64_t pc       = 0;
        int64_t  delta    = 0;
        uint64_t seq      = 0;
        uint64_t prev_seq = 0;
    };

    struct ItEntry {
        uint64_t pc            = 0;
        uint64_t last_cl_addr  = 0;
        uint64_t last_seq      = 0;
        uint16_t frontier_step = 0;
        bool     valid         = false;
    };

    std::array<GhbEntry, GHB_SIZE> ghb_{};
    std::array<ItEntry,  IT_SIZE>  it_{};
    uint64_t                       next_seq_ = 1;

public:
    template <typename EmitFn>
    void Update(uint64_t cl_addr, uint64_t ip, uint64_t /*demand_page*/,
                EmitFn&& /*emit*/)
    {
        ItEntry& it = it_[ip & IT_MASK];
        const bool it_hit = it.valid && it.pc == ip;

        if (!it_hit) {
            it.valid         = true;
            it.pc            = ip;
            it.last_cl_addr  = cl_addr;
            it.last_seq      = 0;
            it.frontier_step = 0;
            return;
        }

        const int64_t delta = static_cast<int64_t>(cl_addr)
                            - static_cast<int64_t>(it.last_cl_addr);
        it.last_cl_addr = cl_addr;

        const uint64_t s = next_seq_++;
        GhbEntry& ge  = ghb_[s % GHB_SIZE];
        ge.pc       = ip;
        ge.delta    = delta;
        ge.seq      = s;
        ge.prev_seq = it.last_seq;
        it.last_seq = s;

        if (it.frontier_step > 0) --it.frontier_step;

        auto window = ReconstructWindow(ip, s);
        (void)window;
    }

private:
    periodicity::PeriodicityBuf<periodicity::WINDOW_SIZE>
    ReconstructWindow(uint64_t pc, uint64_t head_seq) const
    {
        std::array<int64_t, periodicity::WINDOW_SIZE> tmp{};
        std::size_t n = 0;

        uint64_t s = head_seq;
        while (s != 0 && n < periodicity::WINDOW_SIZE) {
            const GhbEntry& e = ghb_[s % GHB_SIZE];
            if (e.seq != s || e.pc != pc) break;
            tmp[n++] = e.delta;
            s = e.prev_seq;
        }

        periodicity::PeriodicityBuf<periodicity::WINDOW_SIZE> w;
        for (std::size_t i = n; i-- > 0;)
            w.Insert(tmp[i]);
        return w;
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
