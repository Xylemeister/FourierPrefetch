# Fourier Prefetcher — Improvement Plan

> **Prerequisite**: Fix all bugs documented in `REVIEW.md` first (broken `DominantBin`,
> disabled `SpectralPurity` gate, N=32 leakage, global tracker). Everything below
> assumes a correct baseline implementation.

---

## Benchmark-Driven Gap Analysis

Comparison against BERTI on 9 ChampSim traces (100M simulation instructions):

| Workload | BERTI IPC | Fourier IPC | Gap |
|---|---|---|---|
| 602.gcc_s-734B | 3.315 | 2.840 | −14% |
| 602.gcc_s-1850B | 3.471 | 2.686 | −23% |
| 602.gcc_s-2226B | 1.232 | 1.241 | ≈ tie |
| 602.gcc_s-2375B | 2.408 | 2.366 | −2% |
| 649.fotonik3d_s-1B | 3.393 | 3.393 | tie |
| 649.fotonik3d_s-1176B | 3.028 | 1.949 | **−36%** |
| 649.fotonik3d_s-7084B | 1.917 | 1.664 | −13% |
| 649.fotonik3d_s-8225B | 3.002 | 1.948 | **−35%** |
| 649.fotonik3d_s-10881B | 1.464 | 1.235 | −16% |

**Useful prefetch counts reveal three distinct failure modes:**

| Workload | BERTI issued | BERTI useful | Fourier issued | Fourier useful |
|---|---|---|---|---|
| gcc_s-1850B | 12,162,745 | 1,554,426 | **11,466** | 326 |
| gcc_s-2226B | 12,655,649 | 735,579 | **6,672** | 0 |
| fotonik3d-10881B | 3,248,572 | 1,004,859 | **0** | 0 |
| fotonik3d-1176B | 5,358,137 | 1,229,625 (23%) | 3,569,635 | 116,222 (**3%**) |
| fotonik3d-8225B | 5,279,337 | 1,212,397 (23%) | 3,567,704 | 114,859 (**3%**) |

**Three root causes, in order of impact:**

1. **Coverage failure** — Fourier issues near-zero prefetches on gcc_s traces where
   patterns are non-stationary. DFT requires the access stream to be periodic across
   the entire `WINDOW_SIZE=64` window; when the stride shifts mid-window the spectrum
   goes flat and `SpectralPurity` drops below threshold → nothing is issued.

2. **No timeliness** — `getSingleDelta(period)` in `transform.h:134` issues exactly
   one prefetch, one period ahead, regardless of memory latency. For fotonik3d-1176B
   the L2C miss latency is **220 cycles** but the stride interval is much shorter —
   the prefetch arrives too late to hide the miss. BERTI explicitly measures latency
   and scales its lookahead depth accordingly.

3. **Precision failure** — For fotonik3d-1176B and 8225B, Fourier issues ~3.5M
   prefetches but only 3% are useful (vs BERTI's 23%). Spurious DFT harmonics at
   `2f, 3f, ...` generate wrong addresses alongside the real stride.

---

## Priority 1 — Close Coverage Failures

### A. Sliding-window STFT (non-stationary patterns)

The current design applies DFT to the entire `WINDOW_SIZE` history. When the program
phase changes (common in gcc_s which is a compiler with many short loops), the old
deltas corrupt the spectrum. Replace the fixed window with an **exponential decay**
applied to the running DFT bins already computed by `InsertTransform`:

```cpp
// In InsertTransform, after updating re_[k] / im_[k]:
re_[k] *= DECAY;   // e.g. DECAY = 0.97
im_[k] *= DECAY;
bin_[k] = re_[k] * re_[k] + im_[k] * im_[k];
```

With `DECAY = 0.97` a sample 64 accesses old has weight `0.97^64 ≈ 0.14` — recent
accesses dominate. This gives the effect of a short-time Fourier transform without
any window management overhead. The `InsertTransform` recurrence in `transform.h:79`
already supports this — just multiply after the twiddle rotation.

Tune `DECAY` as the key sensitivity knob: lower = faster adaptation, higher = more
noise rejection.

### B. Flat-spectrum fallback

When `SpectralPurity()` is below `SPECTRAL_THRESHOLD` for a mature entry, currently
nothing is issued. For gcc_s-1850B this means 12M BERTI prefetches vs near-zero from
Fourier. Add a **stride fallback** at the point in `FourierPrefetchV1::update` where
the purity check fails:

```cpp
if (found->buf.SpectralPurity() > SPECTRAL_THRESHOLD) {
    // ... existing Fourier path ...
} else {
    // fallback: check if last 4 deltas are identical
    auto deltas = found->buf.GetCyclicDeltas(4);
    if (deltas[0] == deltas[1] && deltas[1] == deltas[2] && deltas[2] == deltas[3]
        && deltas[0] != 0) {
        pending_delta_ = deltas[0];
    }
}
```

This adds zero hardware cost (the buffer already exists) and eliminates the
hard-zero cases on non-periodic but strided workloads.

### C. Harmonic suppression

The DFT of a stride-K pattern has peaks at `f, 2f, 3f, ...`. Without suppression,
`DominantBin()` can return a harmonic instead of the fundamental, causing the period
estimate to be halved and the prefetch address to be wrong. After finding the dominant
bin `k*`, zero out harmonic bins before returning:

```cpp
std::size_t DominantBin() const {
    std::size_t best = MinBin();
    for (std::size_t k = MinBin() + 1; k < EFFECTIVE_SIZE; ++k)
        if (bin_[k] > bin_[best]) best = k;
    return best;
}

// call this before SpectralPurity() to suppress harmonics
void SuppressHarmonics(std::size_t fundamental) {
    for (std::size_t k = 2 * fundamental; k < EFFECTIVE_SIZE; k += fundamental)
        bin_[k] = 0.0;
}
```

The correct call sequence in `FourierPrefetchV1::update` becomes:
```cpp
found->buf.SuppressHarmonics(found->buf.DominantBin());
if (found->buf.SpectralPurity() > SPECTRAL_THRESHOLD) { ... }
```

This is why fotonik3d precision is 3% not 23% — fix this before anything else for
those traces.

---

## Priority 2 — Timeliness

This is the single largest remaining gap for fotonik3d after fixing coverage.

### D. Per-PC latency tracking and depth scaling

Add a latency EWMA to `tracker_entry`. Use `prefetcher_cache_fill` (currently a
no-op in `fourier.cc:32`) to close the loop:

```cpp
struct tracker_entry {
    // ... existing fields ...
    uint64_t miss_cycle  = 0;   // cycle when the demand miss was issued
    double   avg_lat     = 80.0; // EWMA of observed fill latency, initialised to ~L2
};
```

In `prefetcher_cache_operate`, on a cache miss (`cache_hit == 0`), record
`entry.miss_cycle = current_cycle`.

In `prefetcher_cache_fill`, look up the IP from metadata and compute:
```cpp
double lat = static_cast<double>(current_cycle - entry.miss_cycle);
entry.avg_lat = 0.875 * entry.avg_lat + 0.125 * lat;
```

Then in the prefetch issue path, scale the number of prefetch depths emitted:
```cpp
// stride_interval = average cycles between accesses to this IP
// depth = how many periods ahead we need to prefetch
std::size_t depth = std::max(1UL,
    static_cast<std::size_t>(std::ceil(entry.avg_lat / stride_interval)));
depth = std::min(depth, MAX_PREFETCH_DEPTH);  // cap e.g. at 6

for (std::size_t d = 1; d <= depth; ++d) {
    int64_t pf_delta = d * getSingleDelta(period);
    // page boundary check, then prefetch_line(...)
}
```

For fotonik3d-1176B with avg_lat ≈ 220 cycles this will naturally push Fourier to
issue 3–4 prefetches ahead, matching what BERTI does via its explicit latency table.

The `stride_interval` can be approximated from the access count and cycle count
tracked per entry, or estimated from `avg_lat / PREFETCH_DEGREE` as a starting point.

---

## Priority 3 — Beyond BERTI's Capabilities

These improvements cannot be replicated by BERTI's delta-correlation design.

### E. Phase-aware issue timing

`InsertTransform` in `transform.h:79` maintains `re_[k]` and `im_[k]` — the full
complex DFT output. The phase of the dominant bin `k*` is:

```cpp
double phase = std::atan2(im_[k_star], re_[k_star]);
```

This tells you *where in the current period* you are, not just what the period is.
Use this to predict the exact access-count offset until the next cycle boundary,
and delay issuing the prefetch until `(lat - cycles_until_next_access)` cycles remain.

Concretely: instead of always issuing immediately on the current access, compute
`cycles_until_boundary = (2π - current_phase) / (2π / period) * stride_interval_cycles`
and store the prefetch in a **delayed issue queue** (indexed by target cycle). This
makes the prefetch arrive right as the next demand access would have been a miss.

No existing prefetcher exploits the imaginary part of the DFT for timing — this is a
genuine novelty.

### F. Multi-component prefetching via top-K reconstruction

Currently `getSingleDelta(period)` returns one stride from the dominant frequency.
Real access streams often have two overlapping strides simultaneously (e.g. a loop
over a struct with two fields at offsets 0 and 64). The DFT captures both as separate
bins. Use the top-2 bins to reconstruct two independent delta sequences and prefetch
both:

```cpp
struct SpectralPeak { std::size_t bin; double power; };

std::array<SpectralPeak, 2> TopTwoPeaks() const {
    // find top-2 non-harmonic bins
}
```

For each peak above a minimum power threshold, issue a separate prefetch stream.
BERTI can only track one delta at a time per PC; Fourier with top-K is structurally
superior for overlapping strides.

### G. Spectral entropy as a dynamic confidence gate

Replace the fixed `SPECTRAL_THRESHOLD` with **spectral entropy**:

```cpp
double SpectralEntropy() const {
    double total = 0.0;
    for (auto p : bin_) total += p;
    if (total == 0.0) return 1.0;
    double H = 0.0;
    for (auto p : bin_) {
        if (p > 0.0) { double r = p / total; H -= r * std::log2(r); }
    }
    return H / std::log2(static_cast<double>(EFFECTIVE_SIZE)); // normalised [0,1]
}
```

Low entropy (< 0.3) → concentrated spectrum → high confidence → prefetch
aggressively (large depth, issue to L1).  
High entropy (> 0.7) → flat spectrum → low confidence → reduce depth, fall back to
stride, or skip entirely.

This is a single, principled knob that unifies coverage and accuracy tradeoffs
adaptively per-PC and per-phase, without hard-coded thresholds that need retuning per
workload.

---

## Idea Space (longer-term)

### Idea 1 — SMS-Fourier hybrid (spatial, page-based)

Existing spatial prefetchers (SMS, DSPatch) record which cache-line offsets within a
4 KB page are accessed together, then replay the pattern on the next trigger. The
pattern is stored as a 64-bit bitvector (one bit per cache line in the page).

Replace the bitvector with a **frequency spectrum over page offsets**: given a
sequence of cache-line offsets visited within a page (0–63), treat those as a
discrete signal and find dominant inter-offset distances with DFT. On a trigger, emit
prefetches for all offsets that fall at multiples of the dominant period from the
trigger offset.

This would handle cases where SMS fails because two pages have the same pattern but
different starting offsets — the Fourier representation is shift-invariant.

### Idea 2 — Fourier-guided SMS pattern pruning

SMS accumulates a 64-bit access bitvector per (region, trigger PC) pair. As the
number of distinct patterns grows, storage pressure forces eviction of useful entries.
Use Fourier analysis to **compress** the bitvector: if the set bits are well-described
by a periodic pattern, store just `(offset, period, count)` instead of the full 64
bits. This reduces the per-entry storage cost by ~4× and allows the SMS table to hold
more distinct program regions simultaneously.

The DFT serves as a lossless compressor for periodic patterns and a lossy compressor
for near-periodic ones — compatible with DSPatch's existing confidence filtering.

---

## Implementation Roadmap

```
Step 1  Fix REVIEW.md bugs                    ← correctness prerequisite
Step 2  Add harmonic suppression (§C)         ← fixes fotonik3d precision 3%→~20%
Step 3  Add flat-spectrum stride fallback (§B) ← fixes gcc_s coverage
Step 4  Add exponential decay STFT (§A)       ← fixes gcc_s non-stationarity
Step 5  Add latency EWMA + depth scaling (§D) ← closes fotonik3d IPC gap
Step 6  Add spectral entropy gating (§G)      ← replace hard threshold
Step 7  Add top-K multi-component (§F)        ← beats Berti on overlapping strides
Step 8  Phase-aware issue timing (§E)         ← novel, benchmark carefully
Step 9  SMS-Fourier (Idea 1/2)                ← full redesign, separate branch
```

Steps 1–5 should be sufficient to match BERTI on the current trace set. Steps 6–8
are where Fourier structurally outperforms BERTI. Steps 1–8 all operate within the
existing `transform.h` + `fourier.cc` architecture.
