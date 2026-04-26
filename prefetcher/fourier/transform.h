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
#include <optional>
#include <vector>
#include "msl/lru_table.h"
#include <sys/types.h>

namespace Transform {

/*****************************************************************************
 *                          Parameters                                       *
 *****************************************************************************/

constexpr std::size_t WINDOW_SIZE        = 64;   // delta history depth / DFT length
constexpr std::size_t TRACKER_SETS       = 256;
constexpr std::size_t TRACKER_WAYS       = 4;
constexpr std::size_t PREFETCH_DEGREE    = 3;  // replaced by full period — MSHR break acts as natural limiter
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


    // Incremental O(K) sliding DFT — exact from the very first sample.
    // buf_ is zero-initialised so outgoing is 0 until a slot is overwritten,
    // equivalent to a zero-padded DFT of the samples seen so far.
    //
    // Recurrence per bin k:
    //   X_k[n] = e^{j·2πk/N} · (X_k[n-1] + x_new − x_old)
    void InsertTransform(int64_t value)
    {
        double outgoing = static_cast<double>(buf_[head_]);
        Insert(value);

        double d = static_cast<double>(value) - outgoing;
        for (std::size_t k = 0; k < EFFECTIVE_SIZE; ++k) {
            double r = re_[k] + d;
            double i = im_[k];
            double c = CosTable_[k][1];   // cos(2πk/N) — twiddle factor e^{j2πk/N}
            double s = SineTable_[k][1];  // sin(2πk/N)
            re_[k]  = r * c - i * s;
            im_[k]  = r * s + i * c;
            bin_[k] = re_[k] * re_[k] + im_[k] * im_[k];
        }
    }

    bool IsMature() const { return count_ >= SIZE; }

    // Lowest bin whose period fits within the samples seen so far.
    // Bins below this correspond to periods > count_, which are meaningless.
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

    static constexpr std::size_t PeriodOfBin(std::size_t k)
    {
        return (k == 0) ? SIZE : SIZE / k;
    }

    // Fraction of total AC power held by the dominant bin.
    // Range [0, 1]: 1.0 = perfectly periodic, 0.0 = no signal.
    // to check if its just by chance, or there is a strong correlation
    // I guess we are skipping the DC component here
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
        // std::size_t n = std::min(period, count_);
        std::size_t n = period;
        std::vector<int64_t> out(n);
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t pos = (head_ + SIZE - n + i) % SIZE;
            out[i] = buf_[pos];
        }
        return out;
    }

    int64_t getSingleDelta(std::size_t period) const {
        std::size_t pos = (head_ + SIZE - period) % SIZE;
        int64_t prefetch_addr = 0;

        // while (pos != head_){
        //     prefetch_addr += buf_[pos];
        //     pos = (pos + 1) % SIZE;
        // }

        for (size_t i = 0; i < period; i++){
            prefetch_addr += buf_[(pos + i) % SIZE];
        }

        return prefetch_addr;
    }

    double getMean() const {
        double res = 0;

        for (auto& elem: buf_){
            res += elem;
        }

        return res /SIZE;
    }

    // Read-only views
    const std::array<double,  EFFECTIVE_SIZE>& viewBin() const { return bin_; }
    const std::array<int64_t, SIZE>&           viewBuf() const { return buf_; }

private:
    std::size_t head_  = 0;
    std::size_t count_ = 0;

    std::array<int64_t, SIZE>           buf_{};
    std::array<double,  EFFECTIVE_SIZE> bin_{};
    std::array<double,  EFFECTIVE_SIZE> re_{};   // running complex DFT bins (real part)
    std::array<double,  EFFECTIVE_SIZE> im_{};   // running complex DFT bins (imag part)

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

struct PrefetchCandidate {
    int64_t delta;    // cache-line delta from the current address
    bool    fill_l1;  // true → fill to L1, false → fill to L2
};

/*
 * Per-IP tracker table.  On each access:
 *   1. Compute the delta from the previous access for this IP.
 *   2. Insert the delta and recompute the spectrum.
 *   3. If the buffer is mature and SpectralPurity ≥ SPECTRAL_THRESHOLD,
 *      derive the dominant period P and return the last P deltas as
 *      prefetch candidates (up to the full detected period).
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
    std::optional<int64_t> pending_delta_;

public:
    void update(uint64_t cl_addr, uint64_t ip)
    {
        pending_delta_ = std::nullopt;
        auto found = table_.check_hit({ip, cl_addr, {}});

        if (found.has_value()) {
            auto stride = static_cast<int64_t>(cl_addr) - static_cast<int64_t>(found->last_cl_addr);
            // found->buf.Insert(stride);
            // found->buf.Transform();
            found->buf.InsertTransform(stride);
            found->last_cl_addr = cl_addr;

            // if (found->buf.IsMature() && found->buf.SpectralPurity() > SPECTRAL_THRESHOLD) {
            if (found->buf.SpectralPurity() > SPECTRAL_THRESHOLD) {
                auto period = found->buf.PeriodOfBin(found->buf.DominantBin());
                pending_delta_ = found->buf.getSingleDelta(period);
            }

            table_.fill(*found);
        } else {
            table_.fill({ip, cl_addr, {}});
        }
    }

    std::optional<int64_t> issue(uint64_t /*ip*/)
    {
        if constexpr (champsim::debug_print) {
            if (pending_delta_.has_value()) {
                    std::cout << "[ISSUED OFFSET] " << pending_delta_.value() << std::endl;
            }
        }
        return std::exchange(pending_delta_, std::nullopt);
    }
};

} // namespace Transform

#endif // TRANSFORM_H
