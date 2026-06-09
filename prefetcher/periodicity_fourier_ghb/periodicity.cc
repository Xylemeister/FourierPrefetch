#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
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

constexpr double MSHR_DEMOTE_RATIO  = 0.5;   // high-prio: fill at L2 below this
constexpr double MSHR_LOWPRIO_RATIO = 0.25;  // low-prio (L2-only): skip at/above this

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

struct Stats {
    uint64_t updates              = 0;
    uint64_t pc_seed              = 0;
    uint64_t not_mature           = 0;
    uint64_t frontier_saturated   = 0;

    uint64_t predict_called       = 0;
    uint64_t no_signal            = 0;
    uint64_t emitted              = 0;
    uint64_t pred_l2_only         = 0;

    std::array<uint64_t,
        periodicity::PeriodicityBuf<periodicity::WINDOW_SIZE>::PMAX + 1>
        period_hist{};
    std::array<uint64_t, 11> score_hist{};

    uint64_t emit_attempts        = 0;
    uint64_t emit_accepted        = 0;
    uint64_t emit_page_filtered   = 0;
    uint64_t emit_pf_filtered     = 0;

    uint64_t l2_attempts          = 0;
    uint64_t l2_accepted          = 0;

    uint64_t c_total              = 0;
    uint64_t c_useful             = 0;
    uint64_t evicted_unused       = 0;
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
    Stats          stats_{};

public:
    void OnUseful(uint64_t cl_addr)
    {
        if (pf_.MarkUseful(cl_addr)) ++stats_.c_useful;
    }

    void OnEvict(uint64_t cl_addr)
    {
        if (pf_.Evict(cl_addr)) ++stats_.evicted_unused;
    }

    void PrintStats(const std::string& name) const
    {
        std::cout << "\n=== Periodicity-GHB (PC-localized) prefetcher stats (" << name << ") ===\n";
        std::cout << "updates              : " << stats_.updates              << '\n';
        std::cout << "pc_seed (IT miss)    : " << stats_.pc_seed              << '\n';
        std::cout << "not_mature           : " << stats_.not_mature           << '\n';
        std::cout << "frontier_saturated   : " << stats_.frontier_saturated   << '\n';
        std::cout << "predict_called       : " << stats_.predict_called       << '\n';
        std::cout << "no_signal            : " << stats_.no_signal            << '\n';
        std::cout << "emitted              : " << stats_.emitted              << '\n';
        std::cout << "pred_l2_only         : " << stats_.pred_l2_only         << '\n';

        std::cout << "period_hist          :";
        for (std::size_t p = 1; p < stats_.period_hist.size(); ++p)
            std::cout << " p" << p << '=' << stats_.period_hist[p];
        std::cout << '\n';

        std::cout << "score_hist           :";
        for (std::size_t i = 0; i < stats_.score_hist.size(); ++i) {
            const double lo = static_cast<double>(i)     / 10.0;
            const double hi = static_cast<double>(i + 1) / 10.0;
            std::cout << ' ' << lo << '-' << hi << '=' << stats_.score_hist[i];
        }
        std::cout << '\n';

        std::cout << "emit_attempts        : " << stats_.emit_attempts        << '\n';
        std::cout << "emit_accepted        : " << stats_.emit_accepted        << '\n';
        std::cout << "emit_page_filtered   : " << stats_.emit_page_filtered   << '\n';
        std::cout << "emit_pf_filtered     : " << stats_.emit_pf_filtered     << '\n';
        std::cout << "l2_attempts          : " << stats_.l2_attempts          << '\n';
        std::cout << "l2_accepted          : " << stats_.l2_accepted          << '\n';
        std::cout << "c_total              : " << stats_.c_total              << '\n';
        std::cout << "c_useful             : " << stats_.c_useful             << '\n';
        std::cout << "evicted_unused       : " << stats_.evicted_unused       << '\n';
        if (stats_.c_total > 0) {
            const double acc = static_cast<double>(stats_.c_useful)
                             / static_cast<double>(stats_.c_total);
            std::cout << "accuracy (C_useful/C_total) : " << acc << '\n';
        }
    }

    template <typename EmitFn>
    void Update(uint64_t cl_addr, uint64_t ip, uint64_t demand_page,
                EmitFn&& emit)
    {
        ++stats_.updates;

        ItEntry& it = it_[ip & IT_MASK];
        const bool it_hit = it.valid && it.pc == ip;

        if (!it_hit) {
            ++stats_.pc_seed;
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

        if (!window.IsMature()) {
            ++stats_.not_mature;
        } else if (it.frontier_step >= periodicity::MAX_FRONTIER_STEPS) {
            ++stats_.frontier_saturated;
        } else {
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
        ++stats_.predict_called;
        const auto pred = window.Predict(
            periodicity::PERIODICITY_THRESHOLD,
            periodicity::PERIODICITY_L2_THRESHOLD,
            periodicity::PREFETCH_DEPTH,
            frontier_step);

        BumpScoreHist(pred.best_score);

        if (pred.period == 0) {
            ++stats_.no_signal;
            return;
        }
        BumpPeriodHist(pred.period);

        if (pred.l2_only) {
            ++stats_.pred_l2_only;
            for (auto it = pred.deltas.rbegin();
                 it != pred.deltas.rend(); ++it) {
                const int64_t d = *it;
                if (d == 0) continue;

                const uint64_t pf_cl   = static_cast<uint64_t>(
                                             static_cast<int64_t>(cl_addr) + d);
                const uint64_t pf_addr = pf_cl << LOG2_BLOCK_SIZE;

                if ((pf_addr >> LOG2_PAGE_SIZE) != demand_page) {
                    ++stats_.emit_page_filtered;
                    continue;
                }
                if (pf_.Contains(pf_cl)) {
                    ++stats_.emit_pf_filtered;
                    continue;
                }

                ++stats_.l2_attempts;
                if (emit(pf_addr, true)) {
                    ++stats_.l2_accepted;
                } else {
                    break;
                }
            }
            return;
        }

        if (!pred.deltas.empty()) ++stats_.emitted;

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
                ++stats_.emit_page_filtered;
                handled();
                continue;
            }
            if (pf_.Contains(pf_cl)) {
                ++stats_.emit_pf_filtered;
                handled();
                continue;
            }

            ++stats_.emit_attempts;
            if (emit(pf_addr, false)) {
                ++stats_.emit_accepted;
                pf_.Insert(pf_cl);
                ++stats_.c_total;
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

    void BumpPeriodHist(std::size_t period)
    {
        if (period < stats_.period_hist.size())
            ++stats_.period_hist[period];
    }

    void BumpScoreHist(double score)
    {
        long bucket = static_cast<long>(std::floor(score * 10.0));
        if (bucket < 0) bucket = 0;
        if (bucket > 10) bucket = 10;
        ++stats_.score_hist[static_cast<std::size_t>(bucket)];
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
                                         bool useful_prefetch,
                                         uint8_t /*type*/,
                                         uint32_t metadata_in)
{
    const uint64_t cl_addr     = addr >> LOG2_BLOCK_SIZE;
    const uint64_t demand_page = addr >> LOG2_PAGE_SIZE;

    auto& engine = engines[this];

    if (useful_prefetch) engine.OnUseful(cl_addr);

    engine.Update(cl_addr, ip, demand_page,
        [this](uint64_t pf_addr, bool l2_only) -> bool {
            const double mshr = this->get_mshr_occupancy_ratio();
            if (l2_only) {
                // Lower-priority (weak) prefetches get a stricter MSHR gate:
                // skip them once the MSHR is busier than MSHR_LOWPRIO_RATIO so
                // they don't crowd out demand misses or high-prio prefetches.
                if (mshr >= MSHR_LOWPRIO_RATIO) return false;
                return prefetch_line(pf_addr, false, 0);  // always LLC-only
            }
            const bool fill_here = mshr < MSHR_DEMOTE_RATIO;
            return prefetch_line(pf_addr, fill_here, 0);
        });

    return metadata_in;
}

void CACHE::prefetcher_cycle_operate() {}

uint32_t CACHE::prefetcher_cache_fill(uint64_t /*addr*/, uint32_t /*set*/,
                                      uint32_t /*way*/, uint8_t /*prefetch*/,
                                      uint64_t evicted_addr,
                                      uint32_t metadata_in)
{
    if (evicted_addr != 0)
        engines[this].OnEvict(evicted_addr >> LOG2_BLOCK_SIZE);
    return metadata_in;
}

void CACHE::prefetcher_final_stats()
{
    auto it = engines.find(this);
    if (it != engines.end()) it->second.PrintStats(this->NAME);
}
