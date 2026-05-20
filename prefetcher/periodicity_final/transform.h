#ifndef PERIODICITY_FINAL_TRANSFORM_H
#define PERIODICITY_FINAL_TRANSFORM_H

/*
 * Periodicity-transform kernel for the Periodicity Prefetcher.
 *
 *   1. Maintain a ring buffer of up to N cache-line deltas. Prediction
 *      starts once count_ ≥ MIN_CYCLES samples are present — the buffer
 *      need not be full. Closes the per-IP warmup gap that would
 *      otherwise stall coverage for ~N demand accesses on every new IP.
 *   2. For each period p ∈ [1, count_ / MIN_CYCLES], project the active
 *      portion of the window onto P_p (the period-p subspace, Sethares
 *      & Staley 1999): each residue class r ∈ [0,p) is replaced by its
 *      in-class mean.
 *   3. Score s_p = ||π_p(x)||² / ||x||²  ∈ [0, 1].
 *      Adjusted score s_p* = s_p − (p−1)/count_.
 *   4. Pick argmax(s_p*); ties go to the smallest p (since
 *      P_p ⊆ P_{kp}, the smallest tied p is the fundamental).
 *   5. If s_p* clears the threshold, replay the last `period` raw
 *      deltas as cumulative prefetch offsets.
 *
 * Why the (p−1)/count_ adjustment:
 *   P_p has p free parameters (one mean per residue class). Under
 *   i.i.d. random data with zero mean, E[s_p] ≈ p/count_ — a bigger
 *   subspace over-fits noise for free, and the available sample budget
 *   is count_, not N. Anchoring at p=1 (constant model already eats
 *   1/count_), the subtraction centres pure noise near 1/count_ and
 *   keeps a real period-p signal near 1 − (p−1)/count_, which still
 *   clears the threshold once a modest number of cycles is in the
 *   window.
 *
 * Why per-p maturity instead of "buffer full":
 *   A clean period-p hypothesis needs only MIN_CYCLES · p samples to
 *   be statistically meaningful. Requiring count_ = N stalls a stride
 *   (p=1) detection by N demand accesses for no mathematical reason.
 *   The p-scan is capped at count_ / MIN_CYCLES so we never score a
 *   period we can't see two cycles of, and the (p−1)/count_ penalty
 *   keeps small windows from being seduced by large p.
 *
 * Why no DC special case:
 *   A constant-stride stream is just a period-1 signal. P_1 is the
 *   1-dim subspace of constants, so a pure stride gives s_1 = 1.0
 *   with zero DoF penalty, and `period = 1` cumulative replay
 *   produces the stride-1, stride-2, stride-3, stride-4 lookahead.
 *   Folding DC into the period scan removes a special case, a
 *   separate emission policy, and an entire bug class (DC stealing
 *   periodic signals via the mean).
 *
 * Cost: O(count_ · p_max) per query, O(1) per Insert.
 */

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace periodicity {

constexpr std::size_t WINDOW_SIZE    = 16;
constexpr std::size_t PREFETCH_DEPTH = 4;

// Soft cap on per-IP lookahead distance, measured in cycle positions
// past the current demand. PREFETCH_DEPTH lines are emitted per demand;
// each demand consumes one cycle position. With cap = 4·PREFETCH_DEPTH
// the steady-state lookahead horizon is ~PREFETCH_DEPTH·(1 - 1/depth)^-1
// cache-line positions ahead of demand before saturation halts new
// emission, matching stride-style frontier behaviour.
constexpr std::size_t MAX_FRONTIER_STEPS = 4 * PREFETCH_DEPTH;

// A period p is considered if at least MIN_CYCLES full cycles fit
// in the window. With N=16 and MIN_CYCLES=2, MAX_PERIOD = 8.
constexpr std::size_t MIN_CYCLES = 2;
constexpr std::size_t MAX_PERIOD = WINDOW_SIZE / MIN_CYCLES;

// Acceptance gates on the DoF-adjusted score. A clean period-p signal
// scores 1 − (p−1)/N here; at PMAX=8 that's 0.56, so the L1 gate of
// 0.50 admits any clean integer period in [1, PMAX].
//
// PERIODICITY_L2_THRESHOLD opens a marginal-signal band [L2, L1) whose
// best period is emitted with fill_here=false → the line fills L2 only.
// Trade-off: L1 capacity is untouched (so wrong guesses don't pollute
// the hot path), and a partially-periodic IP still hides ~80 cycles of
// DRAM latency when we guess right. Set L2 above the random-noise
// floor of ~1/count (~0.06 at N=16) so noise rarely clears it.
constexpr double PERIODICITY_THRESHOLD    = 0.5;
constexpr double PERIODICITY_L2_THRESHOLD = 0.3;

// Per-call cap on non-primary candidates whose deltas are unioned with
// the primary's. The top (1 + MAX_SECONDARY_CANDIDATES) periods by
// adjusted score are considered; secondaries are kept only if they
// clear PERIODICITY_THRESHOLD on their own. Used to widen coverage on
// IPs whose access pattern carries multiple coexisting periods (e.g.
// bwaves, where p=4 and p=5 commonly co-fire on the same IP).
constexpr std::size_t MAX_SECONDARY_CANDIDATES = 1;

struct Prediction {
    std::vector<int64_t> deltas;        // primary period's cumulative offsets
    std::vector<int64_t> extra_deltas;  // secondary candidates' offsets (dedup'd)
    std::size_t          period      = 0;     // 0 → no usable period
    double               best_score  = 0.0;   // max adjusted score
    std::size_t          candidates  = 0;     // # periods clearing L1 threshold
    bool                 l2_only     = false; // marginal: route primary → L2
};

template <std::size_t N>
class PeriodicityBuf
{
public:
    static constexpr std::size_t PMAX = N / MIN_CYCLES;

    void Insert(int64_t value)
    {
        buf_[head_] = value;
        head_       = (head_ + 1) % N;
        if (count_ < N) ++count_;
    }

    // Lowest count_ at which any period (p=1) is statistically meaningful.
    // Predict() further caps p at count_ / MIN_CYCLES inside its scan.
    bool IsMature() const { return count_ >= MIN_CYCLES; }

    // Scan p ∈ [1, PMAX]; populate `deltas` (primary), `extra_deltas`
    // (top-K secondaries that also clear the L1 threshold), and the
    // l2_only flag (best score in the marginal band [L2, L1)).
    //
    // start_step shifts the primary replay so the cumulative offsets
    // land `start_step` cycle positions beyond `cl_addr` instead of
    // starting at the immediate next line. Secondaries always start at
    // step 0 — their cycle positions are in a different period than
    // the frontier-tracked primary, so they're best-effort and the
    // recent-ring filter handles any address overlap.
    Prediction Predict(double      threshold     = PERIODICITY_THRESHOLD,
                       double      l2_threshold  = PERIODICITY_L2_THRESHOLD,
                       std::size_t depth         = PREFETCH_DEPTH,
                       std::size_t start_step    = 0,
                       std::size_t extra_max     = MAX_SECONDARY_CANDIDATES) const
    {
        Prediction out;
        if (depth == 0 || !IsMature()) return out;
        const double total = TotalEnergy();
        if (total <= 0.0) return out;

        struct Scored { double s_adj; std::size_t p; };
        std::array<Scored, PMAX> scored{};
        std::size_t n = 0;

        // Cap p at count_/MIN_CYCLES so we never score a period we
        // haven't observed MIN_CYCLES full cycles of.
        const std::size_t p_max =
            std::min<std::size_t>(PMAX, count_ / MIN_CYCLES);
        for (std::size_t p = 1; p <= p_max; ++p) {
            const double s_raw = ProjectionEnergy(p) / total;
            const double s_adj = s_raw - static_cast<double>(p - 1)
                                       / static_cast<double>(count_);
            if (s_adj >= threshold) ++out.candidates;
            scored[n++] = {s_adj, p};
        }
        if (n == 0) return out;

        // Order the top (1 + extra_max) by adjusted score; ties go to
        // the smaller p (the fundamental, since P_p ⊆ P_{kp}).
        const std::size_t want = std::min<std::size_t>(n, 1 + extra_max);
        std::partial_sort(scored.begin(), scored.begin() + want,
                          scored.begin() + n,
                          [](const Scored& a, const Scored& b) {
                              return a.s_adj > b.s_adj
                                  || (a.s_adj == b.s_adj && a.p < b.p);
                          });

        out.best_score = scored[0].s_adj;
        if (scored[0].s_adj < l2_threshold) return out;

        out.period  = scored[0].p;
        out.deltas  = CyclicReplay(scored[0].p, depth, start_step);
        out.l2_only = (scored[0].s_adj < threshold);
        if (out.l2_only) return out;  // marginal-signal: skip secondaries

        // Append cumulative offsets from each remaining candidate that
        // still clears the L1 threshold. Dedup against `deltas` so we
        // never emit the same cycle-position twice in one call; the
        // recent-ring filter handles cross-call overlap.
        for (std::size_t k = 1; k < want; ++k) {
            if (scored[k].s_adj < threshold) break;
            if (scored[k].p == scored[0].p) continue;
            const auto sd = CyclicReplay(scored[k].p, depth, /*start_step=*/0);
            for (int64_t d : sd) {
                bool dup = false;
                for (int64_t dp : out.deltas) {
                    if (dp == d) { dup = true; break; }
                }
                if (!dup) out.extra_deltas.push_back(d);
            }
        }

        return out;
    }

private:
    std::array<int64_t, N> buf_{};
    std::size_t            head_  = 0;
    std::size_t            count_ = 0;

    // Chronological accessor over the full ring: At(0) is the slot
    // head_ points at, At(N-1) the slot just before. When count_ = N
    // this is oldest→newest. When count_ < N the first N−count_ slots
    // are zero-initialised; use AtValid for energy computations and
    // reserve At for CyclicReplay, which only reads At(N-period..N-1)
    // — the newest `period` samples — and is correct regardless of
    // count_ as long as count_ ≥ period (guaranteed by the p_max cap).
    int64_t At(std::size_t i) const { return buf_[(head_ + i) % N]; }

    // Chronological accessor restricted to valid samples.
    // j ∈ [0, count_-1]; AtValid(0) is the oldest valid sample,
    // AtValid(count_-1) the newest. Required for energy sums so that
    // zero-initialised unwritten slots don't dilute class means.
    int64_t AtValid(std::size_t j) const
    {
        return buf_[(head_ + N - count_ + j) % N];
    }

    double TotalEnergy() const
    {
        double s = 0.0;
        for (std::size_t j = 0; j < count_; ++j) {
            const double x = static_cast<double>(AtValid(j));
            s += x * x;
        }
        return s;
    }

    // ||π_p(x)||² = Σ_r (sum_r)² / count_r, where sum_r and count_r are
    // the within-class running totals. Closed form: π_p replaces each
    // sample with its class mean, so the L2 norm is Σ_r count_r · mean_r²
    // = Σ_r sum_r² / count_r.
    double ProjectionEnergy(std::size_t p) const
    {
        std::array<double, PMAX>        sums{};
        std::array<std::uint32_t, PMAX> counts{};
        for (std::size_t j = 0; j < count_; ++j) {
            const std::size_t r = j % p;
            sums[r] += static_cast<double>(AtValid(j));
            ++counts[r];
        }
        double e = 0.0;
        for (std::size_t r = 0; r < p; ++r) {
            if (counts[r] == 0) continue;
            e += sums[r] * sums[r] / static_cast<double>(counts[r]);
        }
        return e;
    }

    // Replay the last `period` raw deltas as cumulative offsets. For
    // period=1 this is the constant-stride lookahead (1·μ, 2·μ, …).
    //
    // start_step pre-accumulates the first start_step cycle positions
    // without emitting them; the returned offsets cover cycle positions
    // [start_step, start_step + depth). Used by frontier-tracked
    // emission to skip lines already prefetched while still extending
    // the lookahead horizon by `depth` per call.
    std::vector<int64_t>
    CyclicReplay(std::size_t period, std::size_t depth,
                 std::size_t start_step) const
    {
        std::vector<int64_t> out;
        out.reserve(depth);
        int64_t cumulative = 0;
        for (std::size_t i = 0; i < start_step; ++i) {
            cumulative += At(N - period + (i % period));
        }
        for (std::size_t i = 0; i < depth; ++i) {
            cumulative += At(N - period + ((start_step + i) % period));
            out.push_back(cumulative);
        }
        return out;
    }
};

}  // namespace periodicity

#endif  // PERIODICITY_FINAL_TRANSFORM_H
