#include <array>
#include <cstdint>
#include <unordered_map>

#include "cache.h"
#include "msl/lru_table.h"
#include "transform.h"

namespace {

constexpr std::size_t TRACKER_SETS = 256;
constexpr std::size_t TRACKER_WAYS = 4;

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
    struct Entry {
        uint64_t                                              ip            = 0;
        uint64_t                                              last_cl_addr  = 0;
        periodicity::PeriodicityBuf<periodicity::WINDOW_SIZE> buf{};
        uint16_t                                              frontier_step = 0;

        auto index() const { return ip; }
        auto tag()   const { return ip; }
    };

    champsim::msl::lru_table<Entry> table_{TRACKER_SETS, TRACKER_WAYS};
    PrefetchFilter                  pf_{};

public:
    template <typename EmitFn>
    void Update(uint64_t cl_addr, uint64_t ip, uint64_t demand_page,
                EmitFn&& emit)
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

        if (found->buf.IsMature()
            && found->frontier_step < periodicity::MAX_FRONTIER_STEPS) {
            EmitPredictions(*found, cl_addr, demand_page, emit);
        }

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

    template <typename EmitFn>
    void EmitPredictions(Entry& e, uint64_t cl_addr,
                         uint64_t demand_page, EmitFn& emit)
    {
        const auto pred = e.buf.Predict(
            periodicity::PERIODICITY_THRESHOLD,
            periodicity::PERIODICITY_L2_THRESHOLD,
            periodicity::PREFETCH_DEPTH,
            e.frontier_step);

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

        const std::size_t new_step = static_cast<std::size_t>(e.frontier_step)
                                   + advance;
        e.frontier_step = static_cast<uint16_t>(
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
