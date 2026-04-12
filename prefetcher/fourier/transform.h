#ifndef TRANSFORM_H
#define TRANSFORM_H

/*
 * Fourier Prefetcher
 *
 * Observes per-IP cache-line delta streams, applies a DFT to identify
 * dominant spatial periodicities, and emits prefetch candidates when
 * spectral purity exceeds a confidence threshold.
 *
 * Core idea:
 *   - Maintain a sliding window of WINDOW_SIZE deltas per IP.
 *   - After each access, recompute the DFT power spectrum.
 *   - Identify the dominant non-DC bin k → period P = WINDOW_SIZE / k.
 *   - If power concentration (SpectralPurity) is high enough, replay
 *     the last P deltas as prefetch offsets from the current address.
 */

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>
#include "msl/lru_table.h"

namespace Transform {

/*****************************************************************************
 *                          Parameters                                       *
 *****************************************************************************/

constexpr std::size_t WINDOW_SIZE        = 16;   // delta history depth / DFT length
constexpr std::size_t TRACKER_SETS       = 256;
constexpr std::size_t TRACKER_WAYS       = 4;
constexpr std::size_t PREFETCH_DEGREE    = 4;    // max candidates emitted per access
constexpr double      SPECTRAL_THRESHOLD = 0.50; // dominant bin must hold ≥50% of AC power

/*****************************************************************************
 *                          TransformBuf                                     *
 *****************************************************************************/

/*
 * Circular buffer of SIZE int64_t deltas with an attached DFT engine.
 *
 * Only the non-redundant EFFECTIVE_SIZE = SIZE/2 + 1 frequency bins are
 * computed (exploitation of conjugate symmetry for real-valued inputs).
 * Trig tables are precomputed once per instantiation type and shared
 * across all instances of the same SIZE.
 */
template <std::size_t SIZE>
class TransformBuf
{
public:
    static constexpr std::size_t EFFECTIVE_SIZE = SIZE / 2 + 1;
    static constexpr double      PI             = 3.14159265358979323846;

    using TrigTable = std::array<std::array<double, SIZE>, EFFECTIVE_SIZE>;

    TransformBuf() = default;

    void Insert(int64_t delta)
    {
        buf_[head_] = delta;
        head_       = (head_ + 1) % SIZE;
        if (count_ < SIZE) ++count_;
    }

    // Recompute the power spectrum from the current buffer contents.
    // Should be called after each Insert when a transform is desired.
    void Transform()
    {
        for (std::size_t k = 0; k < EFFECTIVE_SIZE; ++k) {
            double re = 0.0, im = 0.0;
            for (std::size_t n = 0; n < SIZE; ++n) {
                std::size_t idx = (head_ + n) % SIZE;  // oldest→newest order
                re += static_cast<double>(buf_[idx]) * CosTable_[k][n];
                im -= static_cast<double>(buf_[idx]) * SineTable_[k][n];
            }
            bin_[k] = re * re + im * im;  // power (magnitude squared)
        }
    }

    // to check if there is zero padding in the buffer or calc
    bool IsMature() const { return count_ >= SIZE; }

    std::size_t DominantBin() const
    {
        std::size_t best = 1;
        for (std::size_t k = 2; k < EFFECTIVE_SIZE; ++k)
            if (bin_[k] > bin_[best]) best = k;
        return best;
    }

    static constexpr std::size_t PeriodOfBin(std::size_t k)
    {
        return (k == 0) ? SIZE : SIZE / k;
    }

    // Fraction of total AC power held by the dominant bin.
    // Range [0, 1]: 1.0 = perfectly periodic, 0.0 = no signal.
    double SpectralPurity() const
    {
        double total = 0.0;
        for (std::size_t k = 1; k < EFFECTIVE_SIZE; ++k) total += bin_[k];
        if (total == 0.0) return 0.0;
        return bin_[DominantBin()] / total;
    }

    std::vector<int64_t> GetCyclicDeltas(std::size_t period) const
    {
        std::size_t n = std::min(period, count_);
        std::vector<int64_t> out(n);
        // The slot just before head_ is the most-recently written entry.
        // Walk back `n` steps from there to get the oldest in the window.
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t pos = (head_ + SIZE - n + i) % SIZE;
            out[i] = buf_[pos];
        }
        return out;
    }

    // Read-only views
    const std::array<double,  EFFECTIVE_SIZE>& viewBin() const { return bin_; }
    const std::array<int64_t, SIZE>&           viewBuf() const { return buf_; }

private:
    std::size_t head_  = 0;
    std::size_t count_ = 0;

    std::array<int64_t, SIZE>          buf_{};
    std::array<double,  EFFECTIVE_SIZE> bin_{};

    // Precomputed trig tables shared across all instances of this SIZE
    static TrigTable ComputeTable(double (*fn)(double))
    {
        TrigTable t{};
        for (std::size_t k = 0; k < EFFECTIVE_SIZE; ++k)
            for (std::size_t n = 0; n < SIZE; ++n)
                t[k][n] = fn((2.0 * PI * static_cast<double>(k) * static_cast<double>(n))
                             / static_cast<double>(SIZE));
        return t;
    }

    static inline const TrigTable SineTable_ = ComputeTable(std::sin);
    static inline const TrigTable CosTable_  = ComputeTable(std::cos);
};

/*****************************************************************************
 *                          PrefetchCandidate                                *
 *****************************************************************************/

struct PrefetchCandidate {
    int64_t delta;    // cache-line delta from the current address
    bool    fill_l1;  // true → fill to L1, false → fill to L2
};

/*****************************************************************************
 *                          FourierPrefetchV1                                *
 *****************************************************************************/

/*
 * Per-IP tracker table.  On each access:
 *   1. Compute the delta from the previous access for this IP.
 *   2. Insert the delta and recompute the spectrum.
 *   3. If the buffer is mature and SpectralPurity ≥ SPECTRAL_THRESHOLD,
 *      derive the dominant period P and return the last P deltas as
 *      prefetch candidates (capped at PREFETCH_DEGREE).
 */
class FourierPrefetchV1
{
    struct tracker_entry {
        uint64_t ip           = 0;
        uint64_t last_cl_addr = 0;
        TransformBuf<WINDOW_SIZE> buf{};

        // lru_table key interface
        auto index() const { return ip; }
        auto tag()   const { return ip; }
    };

    champsim::msl::lru_table<tracker_entry> table_{TRACKER_SETS, TRACKER_WAYS};

public:
    /*
     * Main entry point — call once per demand access.
     *
     * Parameters:
     *   cl_addr : cache-line address (byte address >> LOG2_BLOCK_SIZE)
     *   ip      : instruction pointer of the load/store
     *
     * Returns a (possibly empty) list of PrefetchCandidates.
     * The caller should walk them, compute the prefetch address as
     *   pf_cl += candidate.delta
     * and call prefetch_line(pf_cl << LOG2_BLOCK_SIZE, candidate.fill_l1, …).
     */
    std::vector<PrefetchCandidate> Operate(uint64_t cl_addr, uint64_t ip)
    {
        TransformBuf<WINDOW_SIZE> tmp{};
        auto found = table_.check_hit({ip, cl_addr, tmp});

        if (!found.has_value()) {
            // First time we see this IP — allocate an entry, no prefetch yet.
            table_.fill({ip, cl_addr, tmp});
            return {};
        }

        int64_t delta = static_cast<int64_t>(cl_addr)
                      - static_cast<int64_t>(found->last_cl_addr);
        found->last_cl_addr = cl_addr;
        found->buf.Insert(delta);
        found->buf.Transform();
        table_.fill(*found);

        // Guard: wait until the window is fully populated
        if (!found->buf.IsMature()) return {};

        // Guard: require a clean dominant frequency
        if (found->buf.SpectralPurity() < SPECTRAL_THRESHOLD) return {};

        std::size_t k      = found->buf.DominantBin();
        std::size_t period = TransformBuf<WINDOW_SIZE>::PeriodOfBin(k);
        std::size_t degree = std::min(period, PREFETCH_DEGREE);

        auto cyclic = found->buf.GetCyclicDeltas(period);

        std::vector<PrefetchCandidate> candidates;
        candidates.reserve(degree);
        for (std::size_t i = 0; i < degree && i < cyclic.size(); ++i)
            candidates.push_back({cyclic[i], /*fill_l1=*/true});

        return candidates;
    }

    /*
     * Debug helper — prints the power spectrum and derived candidates.
     * Use in place of Operate when tracing; identical update semantics.
     */
    void DebugOperate(uint64_t cl_addr, uint64_t ip)
    {
        auto candidates = Operate(cl_addr, ip);

        // Re-read the (now-updated) entry to print spectrum state
        TransformBuf<WINDOW_SIZE> tmp{};
        auto found = table_.check_hit({ip, cl_addr, tmp});
        if (!found.has_value()) return;

        std::size_t k      = found->buf.DominantBin();
        std::size_t period = TransformBuf<WINDOW_SIZE>::PeriodOfBin(k);

        std::cout << "[FOURIER]"
                  << " ip="     << std::hex << ip
                  << " cl="     << cl_addr  << std::dec
                  << " mature=" << found->buf.IsMature()
                  << " purity=" << found->buf.SpectralPurity()
                  << " dom_bin=" << k
                  << " period=" << period
                  << " bins:";
        for (auto v : found->buf.viewBin()) std::cout << ' ' << v;
        std::cout << " | pf_deltas:";
        for (auto& c : candidates) std::cout << ' ' << c.delta;
        std::cout << '\n';
    }
};

} // namespace Transform

#endif // TRANSFORM_H
