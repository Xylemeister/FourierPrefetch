# Fourier Prefetcher — Code Review & Improvement Plan

## Summary

The current implementation in `transform.h` has several correctness bugs and design
divergences from `DESIGN.md` that would make it prefetch incorrectly or not at all
in most cases. The issues range from a completely broken `DominantBin()` to disabled
confidence gating and spectral-leakage problems from the wrong buffer size.

---

## Critical Bugs

### 1. `DominantBin()` does not find the maximum-power bin

The commented-out original was correct. The replacement is wrong:

```cpp
// CURRENT (broken)
std::deque<int> top_k;
for (std::size_t k = 1; k < EFFECTIVE_SIZE; ++k) {
    if (bin_[k] != 0) {
        top_k.push_back(k);
        if (top_k.size() > TOP_K)
            top_k.pop_front();   // keeps the last 3 non-zero *indices*, not powers
    }
}
return top_k.front();            // returns the 3rd-to-last non-zero index
```

For a fully populated spectrum (all bins non-zero, which is normal after a real DFT),
this always returns `EFFECTIVE_SIZE - TOP_K = 14`, i.e. bin 14 regardless of which
bin actually dominates. The period inferred is `32/14 ≈ 2`, so the prefetcher
effectively hard-codes period=2 for every stream.

**Fix**: restore the original argmax:
```cpp
std::size_t DominantBin() const {
    std::size_t best = 1;
    for (std::size_t k = 2; k < EFFECTIVE_SIZE; ++k)
        if (bin_[k] > bin_[best]) best = k;
    return best;
}
```

### 2. `SpectralPurity` gate is commented out

```cpp
// if (found->buf.SpectralPurity() < SPECTRAL_THRESHOLD) return {};
```

With this disabled the prefetcher fires on **every** mature access, regardless of
whether the spectrum shows any real periodicity. Random or stride-1 streams will
generate spurious prefetches at whatever period the broken `DominantBin` returns.
Uncomment this line (after fixing `DominantBin`).

---

## Algorithmic Issues

### 3. `WINDOW_SIZE = 32` causes spectral leakage for non-power-of-2 periods

DFT has zero spectral leakage only when the period divides the window length exactly.
With N=32 the clean periods are {1, 2, 4, 8, 16, 32}. A stream with period 3, 6, or
12 spreads its energy across many bins, making `DominantBin` unreliable and
`SpectralPurity` artificially low.

`DESIGN.md` chose N=24 precisely because 24 = LCM(2,3,4,6,8,12), ensuring each
target period maps to a single bin with no leakage. Either:

- **Option A**: Keep full-DFT approach but change `WINDOW_SIZE` to 24 (or any LCM of
  the target periods).
- **Option B**: Switch to the Goertzel approach from `DESIGN.md` and test only the 6
  candidate periods {2,3,4,6,8,12} — cheaper and leak-free by construction.

### 4. Full O(N²) DFT on every access

`Transform()` runs 17 × 32 = 544 multiply-adds after every insert. For a busy L1D
this fires millions of times per benchmark. DESIGN.md suggests:
- Use Goertzel (6 × 24 = 144 ops) on the candidate periods only.
- Run analysis every 4th access, not every access.

This is purely a performance concern for simulation throughput, not simulation
correctness, but it is meaningful.

### 5. `PeriodOfBin` integer truncation for small k

`SIZE / k` truncates for k that does not divide SIZE. For N=32 and k=3, the computed
period is 32/3 = 10 (integer division), but the true period is ~10.67. Replaying 10
deltas instead of the real period produces a phase drift that accumulates. With the
Goertzel + fixed candidate-period approach this is not an issue because periods are
chosen to divide N exactly.

---

## Design Issues

### 6. No confidence / stability tracking

DESIGN.md describes:
- A **stability counter** that requires the same period to be detected N consecutive
  times before committing to it.
- A **confidence counter** (saturating, 0–7) that tracks whether the last prediction
  was correct, gating prefetch issuance.

Without these, a single mis-detected period causes bad prefetches until the entry is
evicted. These counters are the main mechanism against false positives.

### 7. No stride fast-path

Constant-stride streams (the most common case) require a full DFT before any
prefetch is issued. Adding a simple "last 4 deltas identical" fast-path would:
- Start prefetching after 4 accesses (vs 32).
- Save the DFT cost for the dominant workload pattern.
- Degrade gracefully: if the stride changes the DFT path takes over.

### 8. No page boundary check

The walk in `fourier.cc`:
```cpp
walk_cl += static_cast<uint64_t>(c.delta);
prefetch_line(walk_cl << LOG2_BLOCK_SIZE, ...);
```
can cross a page boundary. Hardware prefetchers traditionally do not cross pages
(different physical page, potential fault). Add:
```cpp
if ((pf_addr >> LOG2_PAGE_SIZE) != (addr >> LOG2_PAGE_SIZE)) continue;
```

### 9. No MSHR-aware fill level

All candidates are filled to L1 unconditionally. Under MSHR pressure this degrades
useful-prefetch ratio. `ip_stride` uses:
```cpp
bool fill_l1 = (get_mshr_occupancy_ratio() < 0.5);
```
Should be replicated here.

### 10. Global tracker — not per-cache instance

```cpp
Transform::FourierPrefetchV1 tracker;  // fourier.cc line 7
```

This single global is shared across all CACHE instances that load this prefetcher.
`ip_stride` and `berti` use `std::map<CACHE*, state>` keyed by the cache pointer
so each cache level gets independent state.

---

## Minor Issues

### 11. `#include <iostream>` in a header

`transform.h` includes `<iostream>` for `DebugOperate`. This forces `<iostream>` into
every translation unit that includes the header. Move debug printing to a `.cc` file
or guard it behind a `#ifdef DEBUG`.

### 12. `table_.fill(*found)` write-back pattern

`check_hit` returns a `std::optional<tracker_entry>` (a copy). Mutations to `*found`
are written back via `table_.fill(*found)`. If `lru_table::fill` allocates a new LRU
slot instead of updating in-place, this could evict a different valid entry and leave
a stale copy of the modified entry. Verify `lru_table` semantics or use a pointer /
reference API if available.

### 13. `GetCyclicDeltas` allocates a heap vector per access

`GetCyclicDeltas(period)` returns `std::vector<int64_t>`. With `period` potentially
up to 32 this is a heap allocation on the hot path. Use a `std::array<int64_t, PREFETCH_DEGREE>`
or at least reserve only `degree` elements and pass `degree` directly so the vector
never exceeds 4 elements.

---

## Proposed Revised Design

Below is a revised architecture that addresses all issues above while staying close to
the original concept.

### Parameters

```cpp
constexpr std::size_t BUF_SIZE           = 24;   // LCM(2,3,4,6,8,12) — no leakage
constexpr std::size_t TRACKER_SETS       = 256;
constexpr std::size_t TRACKER_WAYS       = 4;
constexpr std::size_t PREFETCH_DEGREE    = 4;
constexpr std::size_t ANALYSIS_INTERVAL  = 4;    // run Goertzel every 4th insert
constexpr uint8_t     STABILITY_REQUIRED = 3;    // same period 3x before committing
constexpr uint8_t     CONFIDENCE_INIT    = 4;    // starting confidence on period lock
constexpr uint8_t     CONFIDENCE_MAX     = 7;
constexpr uint8_t     CONFIDENCE_MIN_PF  = 3;    // prefetch only if confidence >= 3
constexpr float       POWER_THRESHOLD    = 0.35f;// fraction of total AC power
```

### Candidate Periods and Goertzel Coefficients

```cpp
// Periods that divide BUF_SIZE=24 exactly
constexpr int   CANDIDATE_PERIODS[] = {2, 3, 4, 6, 8, 12};
constexpr int   NUM_CANDIDATES      = 6;
// coeff[i] = 2 * cos(2*pi * (BUF_SIZE/CANDIDATE_PERIODS[i]) / BUF_SIZE)
constexpr float GOERTZEL_COEFF[]    = {-2.0f, -1.0f, 0.0f, 1.0f, 1.4142f, 1.7321f};
```

### Per-IP Entry

```cpp
struct fourier_entry {
    uint64_t ip             = 0;
    uint64_t last_cl_addr   = 0;
    int16_t  delta_buf[BUF_SIZE]{};
    uint8_t  buf_head       = 0;
    uint8_t  buf_count      = 0;
    uint8_t  access_count   = 0;   // for ANALYSIS_INTERVAL gating

    uint8_t  period         = 0;   // detected period (0 = none)
    uint8_t  phase          = 0;   // current position within period
    uint8_t  stability      = 0;   // consecutive detections of same period
    uint8_t  confidence     = 0;   // saturating counter [0, CONFIDENCE_MAX]
    int16_t  pattern[MAX_PERIOD]{};// extracted replay pattern
    int16_t  last_prediction = 0;
    bool     has_prediction  = false;

    auto index() const { return ip; }
    auto tag()   const { return ip; }
};
```

### `prefetcher_cache_operate` Logic

```
1. Lookup IP; allocate if new, return {} on first access.
2. delta = cl_addr - last_cl_addr; update last_cl_addr.
3. Insert delta into circular buffer.

4. CONFIDENCE UPDATE (if has_prediction):
   actual_delta = delta
   if actual_delta == last_prediction: confidence = min(confidence+1, MAX)
   else: confidence = (confidence >= 2) ? confidence - 2 : 0
   if confidence == 0: period = 0, stability = 0

5. STRIDE FAST-PATH:
   if buf_count >= 4 and last 4 deltas are equal and delta != 0:
     predicted_delta = delta; goto PREFETCH

6. GOERTZEL ANALYSIS (every ANALYSIS_INTERVAL accesses, buf_count >= BUF_SIZE):
   For each candidate period P:
     power[P] = goertzel_power(delta_buf, buf_head, BUF_SIZE, coeff[P])
   Normalise powers; best_P = argmax power[]
   If normalised_power[best_P] >= POWER_THRESHOLD:
     If best_P == period: stability++
     Else: period = best_P, stability = 1, confidence = CONFIDENCE_INIT
     If stability >= STABILITY_REQUIRED:
       copy last `period` deltas into pattern[]; phase = 0
   predicted_delta = (period > 0) ? pattern[phase % period] : 0

7. PREFETCH (if period > 0 and confidence >= CONFIDENCE_MIN_PF):
   page_boundary_check(pf_addr, addr)
   fill_l1 = (get_mshr_occupancy_ratio() < 0.5)
   prefetch_line(pf_addr, fill_l1, metadata_in)
   last_prediction = predicted_delta; has_prediction = true
   phase = (phase + 1) % period
```

### Goertzel Implementation

```cpp
float goertzel_power(const int16_t* buf, int head, int N, float coeff) {
    float s0 = 0.0f, s1 = 0.0f;
    for (int i = 0; i < N; ++i) {
        float s_new = static_cast<float>(buf[(head + i) % N]) + coeff * s0 - s1;
        s1 = s0;
        s0 = s_new;
    }
    return s0*s0 + s1*s1 - coeff*s0*s1;
}
```

---

## Summary Table

| Issue | Severity | Fix |
|-------|----------|-----|
| `DominantBin` broken (returns wrong bin) | **Critical** | Restore argmax |
| SpectralPurity gate commented out | **Critical** | Uncomment |
| N=32 spectral leakage | High | Change BUF_SIZE to 24 |
| Full DFT every access | Medium | Goertzel every 4th access |
| No confidence/stability counters | High | Add per DESIGN.md |
| No stride fast-path | Medium | 4-delta identical check |
| No page boundary check | High | Standard guard |
| No MSHR-aware fill level | Low | `get_mshr_occupancy_ratio()` |
| Global tracker (not per-cache) | High | `std::map<CACHE*, state>` |
| `<iostream>` in header | Low | Move to `.cc` |
| Heap vector per access | Low | Stack array or reserve(DEGREE) |
