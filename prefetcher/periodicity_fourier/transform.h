#ifndef PERIODICITY_FOURIER_TRANSFORM_H
#define PERIODICITY_FOURIER_TRANSFORM_H

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace periodicity {

constexpr std::size_t WINDOW_SIZE    = 16;
constexpr std::size_t PREFETCH_DEPTH = 4;

constexpr std::size_t MAX_FRONTIER_STEPS = PREFETCH_DEPTH;

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

private:
    std::array<int64_t, N> buf_{};
    std::size_t            head_  = 0;
    std::size_t            count_ = 0;

    int64_t At(std::size_t i) const { return buf_[(head_ + i) % N]; }

    int64_t AtValid(std::size_t j) const
    {
        return buf_[(head_ + N - count_ + j) % N];
    }
};

}  // namespace periodicity

#endif  // PERIODICITY_FOURIER_TRANSFORM_H
