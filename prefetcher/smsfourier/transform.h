#ifndef PERIODICITY_FOURIER_TRANSFORM_H
#define PERIODICITY_FOURIER_TRANSFORM_H

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace periodicity {

constexpr std::size_t WINDOW_SIZE    = 16;
constexpr std::size_t PREFETCH_DEPTH = 4;

constexpr std::size_t MAX_FRONTIER_STEPS = 4 * PREFETCH_DEPTH;

constexpr std::size_t MIN_CYCLES = 2;
constexpr std::size_t MAX_PERIOD = WINDOW_SIZE / MIN_CYCLES;

constexpr double PERIODICITY_THRESHOLD    = 0.5;
constexpr double PERIODICITY_L2_THRESHOLD = 0.3;

constexpr std::size_t MAX_SECONDARY_CANDIDATES = 1;

struct Prediction {
    std::vector<int64_t> deltas;
    std::vector<int64_t> extra_deltas;
    std::size_t          period      = 0;
    double               best_score  = 0.0;
    std::size_t          candidates  = 0;
    bool                 l2_only     = false;
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

    bool IsMature() const { return count_ >= MIN_CYCLES; }

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
        if (out.l2_only) return out;

        for (std::size_t k = 1; k < want; ++k) {
            if (scored[k].s_adj < threshold) break;
            if (scored[k].p == scored[0].p) continue;
            const auto sd = CyclicReplay(scored[k].p, depth, 0);
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

    int64_t At(std::size_t i) const { return buf_[(head_ + i) % N]; }

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

    double ProjectionEnergy(std::size_t p) const
    {
        double acc = 0.0;
        for (std::size_t k = 0; k < p; ++k) {
            const double omega = 2.0 * M_PI * static_cast<double>(k)
                               / static_cast<double>(p);
            double re = 0.0, im = 0.0;
            for (std::size_t j = 0; j < count_; ++j) {
                const double x     = static_cast<double>(AtValid(j));
                const double phase = omega * static_cast<double>(j);
                re += x * std::cos(phase);
                im -= x * std::sin(phase);
            }
            acc += re * re + im * im;
        }
        return acc / static_cast<double>(count_);
    }

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

#endif  // PERIODICITY_FOURIER_TRANSFORM_H
