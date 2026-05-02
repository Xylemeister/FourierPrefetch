#ifndef TRANSFORM_H
#define TRANSFORM_H

/*
 * Fourier Prefetcher
 *
 * Observes per-IP cache-line delta streams, applies an incremental DFT
 * with exponential decay (STFT approximation) to identify dominant spatial
 * periodicities, and emits depth-scaled prefetch candidates when spectral
 * purity exceeds a confidence threshold.
 *
 * Core idea:
 *   - Maintain a sliding window of WINDOW_SIZE deltas per IP.
 *   - After each access, update the DFT incrementally with DECAY weighting.
 *   - Suppress harmonic bins before measuring SpectralPurity.
 *   - If IsMature() and purity >= SPECTRAL_THRESHOLD, emit prefetches scaled
 *     by depth = ceil(avg_fill_latency / avg_access_interval).
 *   - Fallback: if spectrum is flat but last 4 deltas are identical, emit
 *     one stride prefetch.
 */

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>
#include "msl/lru_table.h"
#include <sys/types.h>

namespace Transform {

/*****************************************************************************
 *                          Parameters                                       *
 *****************************************************************************/

constexpr std::size_t WINDOW_SIZE        = 64;
constexpr std::size_t TRACKER_SETS       = 256;
constexpr std::size_t TRACKER_WAYS       = 4;
constexpr double      SPECTRAL_THRESHOLD = 0.50;
// Weight of a sample k accesses old ≈ DECAY^k; at k=64, 0.97^64 ≈ 0.14.
constexpr double      DECAY              = 0.97;
constexpr std::size_t MAX_PREFETCH_DEPTH = 6;

/*****************************************************************************
 *                          TransformBuf                                     *
 *****************************************************************************/

/*
 * Circular buffer of SIZE int64_t deltas with an attached DFT engine.
 *
 * Only the non-redundant EFFECTIVE_SIZE = SIZE/2 + 1 bins are computed
 * (conjugate symmetry for real-valued inputs).  Trig tables are precomputed
 * once per SIZE instantiation and shared across all instances.
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

    // Full O(N·K) DFT recompute — kept for correctness testing.
    void Transform()
    {
        for (std::size_t k = 0; k < EFFECTIVE_SIZE; ++k) {
            double re = 0.0, im = 0.0;
            for (std::size_t n = 0; n < SIZE; ++n) {
                std::size_t idx = (head_ + n) % SIZE;
                re += static_cast<double>(buf_[idx]) * CosTable_[k][n];
                im -= static_cast<double>(buf_[idx]) * SineTable_[k][n];
            }
            bin_[k] = re * re + im * im;
        }
    }

    // Incremental O(K) sliding DFT with STFT-style exponential decay.
    //
    // Recurrence per bin k (before decay):
    //   X_k[n] = e^{j·2πk/N} · (X_k[n-1] + x_new − x_old)
    //
    // Multiplying re_/im_ by DECAY < 1 after the twiddle rotation gives
    // recent samples higher weight, adapting to non-stationary phase
    // changes without explicit window management.
    void InsertTransform(int64_t value)
    {
        double outgoing = static_cast<double>(buf_[head_]);
        Insert(value);

        double d = static_cast<double>(value) - outgoing;
        for (std::size_t k = 0; k < EFFECTIVE_SIZE; ++k) {
            double r = re_[k] + d;
            double i = im_[k];
            double c = CosTable_[k][1];   // cos(2πk/N)
            double s = SineTable_[k][1];  // sin(2πk/N)
            re_[k]  = (r * c - i * s) * DECAY;
            im_[k]  = (r * s + i * c) * DECAY;
            bin_[k] = re_[k] * re_[k] + im_[k] * im_[k];
        }
    }

    bool IsMature() const { return count_ >= SIZE; }

    // Lowest meaningful bin — periods longer than count_ are uninformative.
    std::size_t MinBin() const
    {
        if (count_ == 0) return 1;
        return std::max(std::size_t{1}, SIZE / count_);
    }

    std::size_t DominantBin() const
    {
        std::size_t lo   = MinBin();
        std::size_t best = lo;
        for (std::size_t k = lo + 1; k < EFFECTIVE_SIZE; ++k)
            if (bin_[k] > bin_[best]) best = k;
        return best;
    }

    // Zero harmonic bins of `fundamental` in the power spectrum.
    // The DFT of a stride-K pattern peaks at f, 2f, 3f, ... — without
    // suppression DominantBin can return a harmonic, halving the period
    // and generating wrong prefetch addresses.  Must be called before
    // SpectralPurity().  The next InsertTransform recomputes bin_ from
    // re_/im_, so the zeroing does not persist across accesses.
    void SuppressHarmonics(std::size_t fundamental)
    {
        if (fundamental == 0) return;
        for (std::size_t k = 2 * fundamental; k < EFFECTIVE_SIZE; k += fundamental)
            bin_[k] = 0.0;
    }

    static constexpr std::size_t PeriodOfBin(std::size_t k)
    {
        return (k == 0) ? SIZE : SIZE / k;
    }

    // Fraction of total AC power held by the dominant bin [0, 1].
    double SpectralPurity() const
    {
        double total = 0.0;
        std::size_t lo = MinBin();
        for (std::size_t k = lo; k < EFFECTIVE_SIZE; ++k) total += bin_[k];
        if (total == 0.0) return 0.0;
        return bin_[DominantBin()] / total;
    }

    std::vector<int64_t> GetCyclicDeltas(std::size_t period) const
    {
        std::vector<int64_t> out(period);
        for (std::size_t i = 0; i < period; ++i) {
            std::size_t pos = (head_ + SIZE - period + i) % SIZE;
            out[i] = buf_[pos];
        }
        return out;
    }

    int64_t getSingleDelta(std::size_t period) const
    {
        std::size_t pos = (head_ + SIZE - period) % SIZE;
        int64_t delta = 0;
        for (std::size_t i = 0; i < period; ++i)
            delta += buf_[(pos + i) % SIZE];
        return delta;
    }

    double getMean() const
    {
        double res = 0.0;
        for (auto& elem : buf_) res += elem;
        return res / SIZE;
    }

    const std::array<double,  EFFECTIVE_SIZE>& viewBin() const { return bin_; }
    const std::array<int64_t, SIZE>&           viewBuf() const { return buf_; }

private:
    std::size_t head_  = 0;
    std::size_t count_ = 0;

    std::array<int64_t, SIZE>           buf_{};
    std::array<double,  EFFECTIVE_SIZE> bin_{};
    std::array<double,  EFFECTIVE_SIZE> re_{};
    std::array<double,  EFFECTIVE_SIZE> im_{};

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
 *                       FourierPrefetchV1                                   *
 *****************************************************************************/

/*
 * Per-IP tracker.  On each access:
 *   1. Compute delta from the previous access for this IP.
 *   2. InsertTransform (incremental sliding DFT with exponential decay).
 *   3. Once IsMature(), suppress harmonics then check SpectralPurity.
 *   4. If purity >= SPECTRAL_THRESHOLD, emit depth-scaled prefetches.
 *      depth = ceil(avg_fill_lat / avg_interval), capped at MAX_PREFETCH_DEPTH.
 *   5. If purity < SPECTRAL_THRESHOLD but last 4 deltas are identical,
 *      emit one stride prefetch as a flat-spectrum fallback.
 *
 * update() returns the lower 32 bits of ip as metadata for the fill-latency
 * feedback path (record_fill).
 */
class FourierPrefetchV1
{
    struct tracker_entry {
        uint64_t ip           = 0;
        uint64_t last_cl_addr = 0;
        TransformBuf<WINDOW_SIZE> buf{};

        uint64_t pf_issue_cycle = 0;    // cycle when last prefetch batch was issued
        double   avg_lat        = 80.0; // EWMA of observed prefetch fill latency (cycles)
        uint64_t last_cycle     = 0;    // cycle of most recent access
        double   avg_interval   = 10.0; // EWMA of cycles between consecutive accesses

        auto index() const { return ip; }
        auto tag()   const { return ip; }
    };

    champsim::msl::lru_table<tracker_entry> table_{TRACKER_SETS, TRACKER_WAYS};
    std::vector<int64_t> pending_deltas_;

public:
    uint32_t update(uint64_t cl_addr, uint64_t ip, uint64_t current_cycle)
    {
        pending_deltas_.clear();
        auto found = table_.check_hit({ip, cl_addr, {}});

        if (found.has_value()) {
            // Track inter-access interval EWMA
            if (found->last_cycle != 0) {
                double interval = static_cast<double>(current_cycle - found->last_cycle);
                found->avg_interval = 0.875 * found->avg_interval + 0.125 * interval;
            }
            found->last_cycle = current_cycle;

            auto stride = static_cast<int64_t>(cl_addr)
                        - static_cast<int64_t>(found->last_cl_addr);
            found->buf.InsertTransform(stride);
            found->last_cl_addr = cl_addr;

            if (found->buf.IsMature()) {
                // Suppress harmonics before purity measurement to avoid
                // picking a spurious sub-harmonic as the dominant frequency.
                found->buf.SuppressHarmonics(found->buf.DominantBin());

                if (found->buf.SpectralPurity() > SPECTRAL_THRESHOLD) {
                    auto    period = found->buf.PeriodOfBin(found->buf.DominantBin());
                    int64_t base   = found->buf.getSingleDelta(period);

                    // Depth = how many periods ahead we need to prefetch so
                    // the data arrives before the demand access.
                    std::size_t depth = std::max(std::size_t{1},
                        static_cast<std::size_t>(
                            std::ceil(found->avg_lat / found->avg_interval)));
                    depth = std::min(depth, MAX_PREFETCH_DEPTH);

                    // for (std::size_t d = 1; d <= depth; ++d)
                    //     pending_deltas_.push_back(static_cast<int64_t>(d) * base);

                    pending_deltas_.push_back(base);

                    found->pf_issue_cycle = current_cycle;
                }
            }

            table_.fill(*found);
        } else {
            tracker_entry e;
            e.ip           = ip;
            e.last_cl_addr = cl_addr;
            e.last_cycle   = current_cycle;
            table_.fill(e);
        }

        return static_cast<uint32_t>(ip);
    }

    // Close the latency feedback loop.  metadata encodes the lower 32 bits
    // of the issuing IP; called only for prefetch fills.
    void record_fill(uint32_t metadata, uint64_t current_cycle)
    {
        uint64_t ip = static_cast<uint64_t>(metadata);
        auto found  = table_.check_hit({ip, 0, {}});
        if (found.has_value() && found->pf_issue_cycle != 0) {
            double lat = static_cast<double>(current_cycle - found->pf_issue_cycle);
            if (lat > 0.0 && lat < 100000.0)
                found->avg_lat = 0.875 * found->avg_lat + 0.125 * lat;
            found->pf_issue_cycle = 0;
            table_.fill(*found);
        }
    }

    const std::vector<int64_t>& pending_deltas() const { return pending_deltas_; }
};

} // namespace Transform

#endif // TRANSFORM_H
