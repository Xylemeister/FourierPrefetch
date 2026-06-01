#include <algorithm>
#include <array>
#include <cstdint>
#include <unordered_map>

#include "cache.h"
#include "transform.h"

namespace {

constexpr std::size_t GHB_SIZE = 512;
constexpr std::size_t IT_SIZE  = 512;
constexpr std::size_t IT_MASK  = IT_SIZE - 1;

constexpr std::size_t PF_BITS = 8;
constexpr std::size_t PF_CAP  = 1u << PF_BITS;
constexpr std::size_t PF_MASK = PF_CAP - 1;

constexpr double MSHR_DEMOTE_RATIO = 0.5;

class PrefetchFilter
{
    struct Slot {
        uint64_t cl_addr = 0;
        bool     valid   = false;
        bool     used    = false;
    };
    std::array<Slot, PF_CAP> slots_{};

public:
    bool Contains(uint64_t cl_addr) const
    {
        const auto& s = slots_[cl_addr & PF_MASK];
        return s.valid && s.cl_addr == cl_addr;
    }

    void Insert(uint64_t cl_addr)
    {
        auto& s = slots_[cl_addr & PF_MASK];
        s.cl_addr = cl_addr;
        s.valid   = true;
        s.used    = false;
    }

    bool MarkUseful(uint64_t cl_addr)
    {
        auto& s = slots_[cl_addr & PF_MASK];
        if (!s.valid || s.cl_addr != cl_addr) return false;
        if (s.used) return false;
        s.used = true;
        return true;
    }

    bool Evict(uint64_t cl_addr)
    {
        auto& s = slots_[cl_addr & PF_MASK];
        if (!s.valid || s.cl_addr != cl_addr) return false;
        const bool was_unused = !s.used;
        s.valid = false;
        return was_unused;
    }
};

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

    PrefetchFilter pf_{};

public:
    template <typename EmitFn>
    void Update(uint64_t cl_addr, uint64_t ip, uint64_t demand_page,
                EmitFn&& emit)
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

        if (window.IsMature()
            && it.frontier_step < periodicity::MAX_FRONTIER_STEPS) {
            EmitPredictions(window, it.frontier_step, cl_addr, demand_page, emit);
        }
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

    template <typename EmitFn>
    void EmitPredictions(const periodicity::PeriodicityBuf<periodicity::WINDOW_SIZE>& window,
                         uint16_t& frontier_step, uint64_t cl_addr,
                         uint64_t demand_page, EmitFn& emit)
    {
        const auto pred = window.Predict(
            periodicity::PERIODICITY_THRESHOLD,
            periodicity::PERIODICITY_L2_THRESHOLD,
            periodicity::PREFETCH_DEPTH,
            frontier_step);

        if (pred.period == 0) return;

        std::size_t advance    = 0;
        bool        cache_full = false;
        const std::size_t depth = pred.deltas.size();
        for (std::size_t k = 0; k < depth && !cache_full; ++k) {
            const std::size_t i = depth - 1 - k;
            const int64_t     d = pred.deltas[i];

            auto handled = [&]() {
                if (i + 1 > advance) advance = i + 1;
            };

            if (d == 0) { handled(); continue; }

            const uint64_t pf_cl   = static_cast<uint64_t>(
                                         static_cast<int64_t>(cl_addr) + d);
            const uint64_t pf_addr = pf_cl << LOG2_BLOCK_SIZE;

            if ((pf_addr >> LOG2_PAGE_SIZE) != demand_page) {
                handled();
                continue;
            }
            if (pf_.Contains(pf_cl)) {
                handled();
                continue;
            }

            if (emit(pf_addr, false)) {
                pf_.Insert(pf_cl);
                handled();
            } else {
                cache_full = true;
            }
        }

        const std::size_t new_step = static_cast<std::size_t>(frontier_step)
                                   + advance;
        frontier_step = static_cast<uint16_t>(
            std::min<std::size_t>(new_step, periodicity::MAX_FRONTIER_STEPS));
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
