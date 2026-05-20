#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>

#include "cache.h"
#include "msl/lru_table.h"
#include "transform.h"

// =============================================================================
// Periodicity Prefetcher
//
// Per IP, project the recent delta stream onto period-p subspaces
// (Sethares & Staley periodicity transform) and emit prefetches at the
// predicted next addresses. Single emission path — constant strides
// are p=1 and flow through the same cyclic-replay machinery as any
// other period.
//
// Side state:
//   - OriginTable maps in-flight prefetched cl_addr back to the IP
//     that issued it, so cache-side useful/useless events can credit
//     the originator and feed the per-IP throttle.
//   - RecentRing (per Entry) suppresses re-emission on overlapping
//     cyclic replays from successive demands (safety net once frontier
//     tracking is in play — see below — but kept because period flips
//     can still produce overlap).
//   - frontier_step (per Entry) tracks how many cycle positions past
//     `cl_addr` the prefetch frontier sits. Each new demand shifts the
//     window by one delta, so frontier_step is decremented; each
//     successful emit advances it. CyclicReplay starts at
//     `start_step = frontier_step`, so every cycle emits PREFETCH_DEPTH
//     *fresh* lines past the frontier instead of re-issuing lines
//     already in flight. MAX_FRONTIER_STEPS caps the horizon — once the
//     frontier sits that far ahead of demand, emission pauses until the
//     demand stream catches up.
// =============================================================================

namespace {

constexpr std::size_t TRACKER_SETS = 256;
constexpr std::size_t TRACKER_WAYS = 4;

constexpr uint32_t THROTTLE_MIN_OBS      = 64;
constexpr uint32_t THROTTLE_USELESS_MULT = 8;
constexpr uint32_t THROTTLE_DECAY_AT     = 128;

constexpr std::size_t RECENT_RING_CAP = 8;

constexpr std::size_t ORIGIN_BITS = 14;
constexpr std::size_t ORIGIN_CAP  = 1u << ORIGIN_BITS;
constexpr std::size_t ORIGIN_MASK = ORIGIN_CAP - 1;

constexpr double MSHR_DEMOTE_RATIO = 0.5;

// ---------------------------------------------------------------------------
// OriginTable — prefetched cl_addr → originating IP
// ---------------------------------------------------------------------------

class OriginTable
{
    struct Slot {
        uint64_t cl_addr = 0;
        uint64_t ip      = 0;
        bool     valid   = false;
    };
    std::array<Slot, ORIGIN_CAP> slots_{};

public:
    void Record(uint64_t cl_addr, uint64_t ip)
    {
        auto& s = slots_[cl_addr & ORIGIN_MASK];
        s.cl_addr = cl_addr;
        s.ip      = ip;
        s.valid   = true;
    }

    std::optional<uint64_t> Take(uint64_t cl_addr)
    {
        auto& s = slots_[cl_addr & ORIGIN_MASK];
        if (!s.valid || s.cl_addr != cl_addr) return std::nullopt;
        s.valid = false;
        return s.ip;
    }
};

// ---------------------------------------------------------------------------
// RecentRing — per-IP set of recently-prefetched cl_addrs
// ---------------------------------------------------------------------------

class RecentRing
{
    std::array<uint64_t, RECENT_RING_CAP> buf_{};
    uint8_t                               head_ = 0;
    uint8_t                               size_ = 0;

public:
    bool Contains(uint64_t cl_addr) const
    {
        for (uint8_t i = 0; i < size_; ++i)
            if (buf_[i] == cl_addr) return true;
        return false;
    }

    void Push(uint64_t cl_addr)
    {
        buf_[head_] = cl_addr;
        head_       = (head_ + 1) % RECENT_RING_CAP;
        if (size_ < RECENT_RING_CAP) ++size_;
    }
};

// ---------------------------------------------------------------------------
// Stats — diagnostic counters, printed in prefetcher_final_stats
// ---------------------------------------------------------------------------
//
//   updates / miss_table / not_mature / throttle_blocks
//       Where demand accesses go. throttle_blocks dominating updates
//       on a low-coverage trace means the per-IP throttle is gating
//       emission too aggressively.
//
//   predict_called / no_signal / emitted
//       Outcome of Predict(). (predict_called − no_signal − emitted)
//       is the silent middle — a period was found but every replay
//       delta was filtered downstream by page/recent.
//
//   period_hist : indexed by period value in [1, PMAX]. p=1 is the
//       constant-stride bucket; high counts here imply a streaming
//       workload. Non-trivial buckets at p=3,5,7 indicate the kernel
//       is finding integer periods that a power-of-two analysis
//       would have missed.
//
//   score_hist : best DoF-adjusted score in [0,1] bucketed into
//       deciles. Diagnostic for the threshold setting: dense activity
//       below the threshold bucket means strong signals are being
//       rejected; activity at the very top means clean periodicity
//       is the rule, not the exception.
//
//   emit_* : prefetch attempt → cache accept vs. page/recent filter.
//   *_credits / *_lost / origins_recorded
//       Origin-table effectiveness. (useful_credits + useless_credits)
//       / origins_recorded is the attribution hit rate; a low ratio
//       means the throttle counters are noisy.
//
// useless_lost is dominated by non-prefetch evictions and is
// informational only.
struct Stats {
    uint64_t updates              = 0;
    uint64_t miss_table           = 0;
    uint64_t not_mature           = 0;
    uint64_t throttle_blocks      = 0;
    uint64_t frontier_saturated   = 0;

    uint64_t predict_called       = 0;
    uint64_t no_signal            = 0;
    uint64_t emitted              = 0;
    uint64_t pred_l2_only         = 0;     // Predict returned l2_only=true
    uint64_t pred_with_extras     = 0;     // Predict returned non-empty extras

    // Indexed by period value; entry 0 is never written.
    std::array<uint64_t,
        periodicity::PeriodicityBuf<periodicity::WINDOW_SIZE>::PMAX + 1>
        period_hist{};

    // Best DoF-adjusted score bucketed into deciles. Index 10 holds
    // the rare exact-1.0 case (clean signal at p=1, zero penalty).
    std::array<uint64_t, 11> score_hist{};

    uint64_t emit_attempts        = 0;
    uint64_t emit_accepted        = 0;
    uint64_t emit_page_filtered   = 0;
    uint64_t emit_recent_filtered = 0;

    uint64_t l2_attempts          = 0;     // marginal-band L2-only emits
    uint64_t l2_accepted          = 0;
    uint64_t extra_attempts       = 0;     // secondary-candidate emits
    uint64_t extra_accepted       = 0;

    uint64_t origins_recorded     = 0;
    uint64_t useful_credits       = 0;
    uint64_t useful_lost          = 0;
    uint64_t useless_credits      = 0;
    uint64_t useless_lost         = 0;
};

// ---------------------------------------------------------------------------
// PeriodicityPrefetcher — engine, one per cache instance
// ---------------------------------------------------------------------------

class PeriodicityPrefetcher
{
    struct Entry {
        uint64_t                                              ip            = 0;
        uint64_t                                              last_cl_addr  = 0;
        periodicity::PeriodicityBuf<periodicity::WINDOW_SIZE> buf{};
        uint32_t                                              useful        = 0;
        uint32_t                                              useless       = 0;
        RecentRing                                            recent{};
        // Cycle positions the prefetch frontier sits past cl_addr.
        // Decremented per buffer shift (demand consumes one position);
        // advanced per successful emit. Capped at MAX_FRONTIER_STEPS.
        uint16_t                                              frontier_step = 0;

        auto index() const { return ip; }
        auto tag()   const { return ip; }

        bool ThrottleBlocks() const
        {
            const uint32_t obs = useful + useless;
            return obs >= THROTTLE_MIN_OBS
                && useless > THROTTLE_USELESS_MULT * useful;
        }
    };

    champsim::msl::lru_table<Entry> table_{TRACKER_SETS, TRACKER_WAYS};
    OriginTable                     origins_{};
    Stats                           stats_{};

public:
    void OnUseful(uint64_t cl_addr)
    {
        if (auto ip = origins_.Take(cl_addr)) {
            ++stats_.useful_credits;
            Credit(*ip, /*useful=*/true);
        } else {
            ++stats_.useful_lost;
        }
    }

    void OnEvict(uint64_t cl_addr)
    {
        if (auto ip = origins_.Take(cl_addr)) {
            ++stats_.useless_credits;
            Credit(*ip, /*useful=*/false);
        } else {
            ++stats_.useless_lost;
        }
    }

    void PrintStats(const std::string& name) const
    {
        std::cout << "\n=== Periodicity prefetcher stats (" << name << ") ===\n";
        std::cout << "updates              : " << stats_.updates              << '\n';
        std::cout << "miss_table (new IPs) : " << stats_.miss_table           << '\n';
        std::cout << "not_mature           : " << stats_.not_mature           << '\n';
        std::cout << "throttle_blocks      : " << stats_.throttle_blocks      << '\n';
        std::cout << "frontier_saturated   : " << stats_.frontier_saturated   << '\n';
        std::cout << "predict_called       : " << stats_.predict_called       << '\n';
        std::cout << "no_signal            : " << stats_.no_signal            << '\n';
        std::cout << "emitted              : " << stats_.emitted              << '\n';
        std::cout << "pred_l2_only         : " << stats_.pred_l2_only         << '\n';
        std::cout << "pred_with_extras     : " << stats_.pred_with_extras     << '\n';

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
        std::cout << "emit_recent_filtered : " << stats_.emit_recent_filtered << '\n';
        std::cout << "l2_attempts          : " << stats_.l2_attempts          << '\n';
        std::cout << "l2_accepted          : " << stats_.l2_accepted          << '\n';
        std::cout << "extra_attempts       : " << stats_.extra_attempts       << '\n';
        std::cout << "extra_accepted       : " << stats_.extra_accepted       << '\n';
        std::cout << "origins_recorded     : " << stats_.origins_recorded     << '\n';
        std::cout << "useful_credits       : " << stats_.useful_credits       << '\n';
        std::cout << "useful_lost          : " << stats_.useful_lost          << '\n';
        std::cout << "useless_credits      : " << stats_.useless_credits      << '\n';
        std::cout << "useless_lost (noisy) : " << stats_.useless_lost         << '\n';
    }

    template <typename EmitFn>
    void Update(uint64_t cl_addr, uint64_t ip, uint64_t demand_page,
                EmitFn&& emit)
    {
        ++stats_.updates;

        Entry probe;
        probe.ip = ip;
        auto found = table_.check_hit(probe);

        if (!found) {
            ++stats_.miss_table;
            Entry seed;
            seed.ip           = ip;
            seed.last_cl_addr = cl_addr;
            table_.fill(seed);
            return;
        }

        AdvanceStream(*found, cl_addr);
        // Buffer shifted by one delta; the frontier moves one cycle
        // position closer to the demand stream.
        if (found->frontier_step > 0) --found->frontier_step;

        if (!found->buf.IsMature()) {
            ++stats_.not_mature;
        } else if (found->ThrottleBlocks()) {
            ++stats_.throttle_blocks;
        } else if (found->frontier_step >= periodicity::MAX_FRONTIER_STEPS) {
            ++stats_.frontier_saturated;
        } else {
            EmitPredictions(*found, cl_addr, ip, demand_page, emit);
        }

        table_.fill(*found);
    }

private:
    void Credit(uint64_t ip, bool useful)
    {
        Entry probe;
        probe.ip = ip;
        auto found = table_.check_hit(probe);
        if (!found) return;

        if (useful) ++found->useful;
        else        ++found->useless;

        if (found->useful + found->useless >= THROTTLE_DECAY_AT) {
            found->useful  >>= 1;
            found->useless >>= 1;
        }
        table_.fill(*found);
    }

    static void AdvanceStream(Entry& e, uint64_t cl_addr)
    {
        const int64_t stride = static_cast<int64_t>(cl_addr)
                             - static_cast<int64_t>(e.last_cl_addr);
        e.buf.Insert(stride);
        e.last_cl_addr = cl_addr;
    }

    template <typename EmitFn>
    void EmitPredictions(Entry& e, uint64_t cl_addr, uint64_t ip,
                         uint64_t demand_page, EmitFn& emit)
    {
        ++stats_.predict_called;
        const auto pred = e.buf.Predict(
            periodicity::PERIODICITY_THRESHOLD,
            periodicity::PERIODICITY_L2_THRESHOLD,
            periodicity::PREFETCH_DEPTH,
            e.frontier_step);

        BumpScoreHist(pred.best_score);

        if (pred.period == 0) {
            ++stats_.no_signal;
            return;
        }
        BumpPeriodHist(pred.period);

        // Marginal-signal path: emit the primary period's deltas to L2
        // only. No frontier advancement (next high-confidence call
        // should still produce L1 emissions at these positions) and no
        // origin recording (L2 fills don't trigger L1's useful_prefetch
        // hook, so attribution would always look "useless" and skew
        // the throttle).
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
                if (e.recent.Contains(pf_cl)) {
                    ++stats_.emit_recent_filtered;
                    continue;
                }

                ++stats_.l2_attempts;
                if (emit(pf_addr, /*l2_only=*/true)) {
                    ++stats_.l2_accepted;
                    e.recent.Push(pf_cl);
                } else {
                    break;
                }
            }
            return;
        }

        if (!pred.deltas.empty())       ++stats_.emitted;
        if (!pred.extra_deltas.empty()) ++stats_.pred_with_extras;

        // Primary L1 path: furthest-first emission with frontier
        // accounting. Under MSHR pressure, the deepest lookahead is the
        // most valuable line to get in flight — it has the longest
        // latency to amortise, and closer lines have more chance to be
        // re-attempted (or demand-issued) on the next cycle. Advance
        // the frontier to the highest cycle position we managed to
        // handle; gaps below it are skipped.
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
            if (e.recent.Contains(pf_cl)) {
                ++stats_.emit_recent_filtered;
                handled();
                continue;
            }

            ++stats_.emit_attempts;
            if (emit(pf_addr, /*l2_only=*/false)) {
                ++stats_.emit_accepted;
                e.recent.Push(pf_cl);
                origins_.Record(pf_cl, ip);
                ++stats_.origins_recorded;
                handled();
            } else {
                cache_full = true;
            }
        }

        const std::size_t new_step = static_cast<std::size_t>(e.frontier_step)
                                   + advance;
        e.frontier_step = static_cast<uint16_t>(
            std::min<std::size_t>(new_step, periodicity::MAX_FRONTIER_STEPS));

        // Secondary candidates: best-effort, no frontier accounting.
        // These already cleared the L1 threshold on their own period,
        // so they go to L1 as well. Furthest-first within the extras
        // list keeps the same MSHR-resilience principle.
        if (cache_full) return;
        for (auto it = pred.extra_deltas.rbegin();
             it != pred.extra_deltas.rend(); ++it) {
            const int64_t d = *it;
            if (d == 0) continue;

            const uint64_t pf_cl   = static_cast<uint64_t>(
                                         static_cast<int64_t>(cl_addr) + d);
            const uint64_t pf_addr = pf_cl << LOG2_BLOCK_SIZE;

            if ((pf_addr >> LOG2_PAGE_SIZE) != demand_page) {
                ++stats_.emit_page_filtered;
                continue;
            }
            if (e.recent.Contains(pf_cl)) {
                ++stats_.emit_recent_filtered;
                continue;
            }

            ++stats_.extra_attempts;
            if (emit(pf_addr, /*l2_only=*/false)) {
                ++stats_.extra_accepted;
                e.recent.Push(pf_cl);
                origins_.Record(pf_cl, ip);
                ++stats_.origins_recorded;
            } else {
                break;
            }
        }
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

// =============================================================================
// ChampSim hooks
// =============================================================================

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
            const bool fill_here = !l2_only
                && this->get_mshr_occupancy_ratio() < MSHR_DEMOTE_RATIO;
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
