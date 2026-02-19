# Fourier Series Prefetcher for ChampSim

## Context

The goal is to build a novel data prefetcher that uses Fourier analysis to detect periodic patterns in memory access streams. Existing prefetchers (ip_stride, berti, spp) use simple stride detection or delta confidence tables, which can miss complex periodic patterns like alternating strides `[4, 8, 4, 8]`, nested loop patterns `[4, 4, 4, -12]`, or longer repeating sequences. A Fourier-based approach can detect the *period* of such patterns and then replay the pattern cyclically to predict future accesses.

The plan is to build this incrementally — first validate the Fourier analysis itself, then wire it into prefetching.

## Design Overview

**Input**: Per-IP sequence of deltas (differences between consecutive cache-line addresses accessed by the same instruction pointer).

**Analysis**: Use the Goertzel algorithm (efficient single-frequency DFT) to detect dominant periods in the delta sequence.

**Prediction**: Once a period P is detected with sufficient confidence, replay the last P deltas cyclically to predict the next access address.

**Stride fast-path**: If the last 4 deltas are identical, skip Fourier analysis and use simple stride prediction (handles the common case cheaply).

## Data Structures

### Per-IP Entry (`fourier_entry`)
```cpp
struct fourier_entry {
    uint64_t ip;                    // instruction pointer (for lookup)
    uint64_t last_cl_addr;          // last cache-line address (to compute deltas)
    int16_t  delta_buf[BUF_SIZE];   // circular buffer of recent deltas (BUF_SIZE=24)
    uint8_t  buf_head;              // write pointer into delta_buf
    uint8_t  buf_count;             // number of valid entries (up to BUF_SIZE)
    int16_t  pattern[MAX_PERIOD];   // replay buffer (last P deltas, MAX_PERIOD=12)
    uint8_t  period;                // detected period (0=none, 2-12)
    uint8_t  phase;                 // current position within the period
    uint8_t  confidence;            // saturating counter 0-7
    uint8_t  stability;             // consecutive detections of same period
    int16_t  last_prediction;       // delta we predicted last time (for verification)
    bool     has_prediction;        // whether last_prediction is valid
};
// ~78 bytes per entry
```

### IP Tracker Table
- `std::map<uint64_t, fourier_entry*>` keyed by hashed IP, with `std::queue` for FIFO eviction (same pattern as BERTI's `bertit` / `bertit_queue` in `prefetcher/berti/berti.h:147-148`)
- Max 256 entries -> ~20 KB total
- Per-cache-instance state via `std::map<CACHE*, fourier_state>` (same pattern as ip_stride's `std::map<CACHE*, tracker>` at `prefetcher/ip_stride/ip_stride.cc:78`)

### Goertzel Coefficient Table
- 6 compile-time constants for candidate periods {2, 3, 4, 6, 8, 12}
- Each stores `coeff = 2 * cos(2*pi*k/N)` where k = N/P, N=24
- Precomputed values: `{-2.0, -1.0, 0.0, 1.0, 1.4142, 1.7321}` (for periods 2, 3, 4, 6, 8, 12)

## Why N=24 and periods {2,3,4,6,8,12}

24 is the LCM of all candidate periods, so each period divides evenly into the buffer length. This eliminates spectral leakage -- a critical correctness issue with DFT. The candidate periods cover common loop nesting depths in real programs.

## Algorithm

### On `prefetcher_cache_operate(addr, ip, cache_hit, ...)`

```
1. Look up IP in tracker map (hash IP first, same as BERTI: ip ^ (ip >> 2) ^ (ip >> 5))
2. If new IP: allocate entry, evict oldest if table full (FIFO)
3. If returning IP: compute delta = current_cl_addr - entry.last_cl_addr
4. Append delta to circular buffer, update last_cl_addr

5. CONFIDENCE CHECK (before predicting):
   If entry has a pending prediction (has_prediction == true):
     If actual_delta == last_prediction: confidence++ (saturate at 7)
     Else: confidence -= 2 (floor at 0)
     If confidence == 0: reset period and stability

6. STRIDE FAST-PATH:
   If last 4 deltas in buffer are identical and delta != 0:
     -> predicted_delta = delta
     -> Skip to step 8 (issue prefetch)

7. FOURIER ANALYSIS (every 4th access when buf_count >= 24):
   For each candidate period P in {2, 3, 4, 6, 8, 12}:
     Run Goertzel on delta_buf for frequency k = N/P
     Compute power |X(k)|^2
   best_period = P with highest power
   If power > POWER_THRESHOLD:
     If best_period == entry.period: stability++ (already tracking this period)
     Else: entry.period = best_period, stability = 1
     If stability >= STABILITY_REQUIRED (3):
       Extract last P deltas from buffer into pattern[], phase = 0
   predicted_delta = pattern[phase % period]

8. PREFETCH (if period is locked OR stride fast-path, and confidence >= 3):
   pf_addr = addr + predicted_delta * BLOCK_SIZE
   Page boundary check: (pf_addr >> LOG2_PAGE_SIZE) == (addr >> LOG2_PAGE_SIZE)
   Fill level: L1 if confidence >= 5 AND get_mshr_occupancy_ratio() < 0.5, else L2
   prefetch_line(pf_addr, fill_this_level, metadata_in)
   Store predicted_delta in last_prediction, set has_prediction = true
   Advance phase
```

### On `prefetcher_cache_fill(...)` -- no-op initially

### On `prefetcher_cycle_operate()` -- no-op

### On `prefetcher_final_stats()` -- print statistics:
- Number of stride-fast-path predictions vs Fourier-based predictions
- Period distribution histogram (how many IPs locked to each period)
- Average confidence, prediction accuracy

## Goertzel Algorithm

```cpp
// For buffer of N samples, test frequency k = N/P
// coeff = 2 * cos(2 * pi * k / N), precomputed per period
float goertzel_power(const int16_t* buf, int head, int N, float coeff) {
    float s1 = 0, s2 = 0;
    for (int i = 0; i < N; i++) {
        float s0 = buf[(head + i) % N] + coeff * s1 - s2;
        s2 = s1;
        s1 = s0;
    }
    return s1*s1 + s2*s2 - coeff*s1*s2;  // |X(k)|^2
}
```

144 multiply-adds total (24 per period x 6 periods). Negligible vs simulation cost.

## Files to Create

```
prefetcher/fourier/
    fourier.h               -- fourier_entry struct, fourier_state class, Goertzel coefficients
    fourier_parameters.h    -- tunable #defines (BUF_SIZE, MAX_PERIOD, thresholds, table size)
    fourier.cc              -- CACHE:: hook implementations
```

## Files to Modify

- `champsim_config.json` -- set `"L1D": { "prefetcher": "fourier" }`

## Implementation Phases

### Phase 1: Skeleton + Delta Buffer + Stride Fast-Path
Create the three files with all hook functions wired up. Implement:
- IP tracker map with FIFO eviction
- Delta circular buffer
- Stride fast-path (if last 4 deltas match, prefetch next stride)
- No Fourier analysis yet

**Validate**: build and run, compare IPC to ip_stride -- should be similar since stride fast-path covers the same patterns.

### Phase 2: Goertzel Period Detection (debug/logging only)
Add Goertzel computation. On every 4th access when buffer is full:
- Compute power for all 6 candidate periods
- Track best period and stability counter
- Print detected periods and power values in `prefetcher_final_stats()`
- **No prefetching from Fourier yet** -- just observe and log

**Validate**: run on a trace, check that stride-4 accesses show period=1 as dominant (constant signal), alternating patterns show period=2, etc.

### Phase 3: Pattern Replay Prefetching
Enable pattern extraction and replay-based prefetching:
- When stability threshold met, extract pattern and issue prefetches
- Add confidence tracking with prediction verification
- Add MSHR-aware fill level and page boundary checks

**Validate**: benchmark IPC against ip_stride and berti on DPC-3 traces.

### Phase 4: Tuning and Extensions (future)
- Tune POWER_THRESHOLD, STABILITY_REQUIRED, confidence thresholds
- Experiment with larger buffers or more candidate periods
- Consider hybrid with Berti for non-periodic patterns

## Config Change

```json
"L1D": {
    "prefetcher": "fourier",
    "prefetch_as_load": false,
    "virtual_prefetch": false,
    "prefetch_activate": "LOAD,PREFETCH"
}
```

## Verification

1. `./config.sh champsim_config.json && make` -- must compile cleanly
2. `bin/champsim --warmup_instructions 200000000 --simulation_instructions 500000000 <trace.xz>`
3. Phase 1: IPC should be close to ip_stride
4. Phase 2: `prefetcher_final_stats()` output should show reasonable period detections
5. Phase 3: IPC should match or improve over ip_stride on periodic workloads

## Reference Files

| File | What to reference |
|------|-------------------|
| `prefetcher/ip_stride/ip_stride.cc` | Per-IP tracking, page boundary check, MSHR-aware fill, `prefetch_line()` usage |
| `prefetcher/berti/berti.h` | `std::map` + `std::queue` FIFO eviction pattern, namespace organization, per-CPU state via vector |
| `prefetcher/berti/berti_parameters.h` | How to organize tunable `#define` constants in a separate header |
| `prefetcher/no/no.cc` | Minimal template showing all required hook signatures |
| `inc/cache.h` | CACHE class API: `prefetch_line()`, `get_mshr_occupancy_ratio()`, `current_cycle`, `NUM_SET`, `NUM_WAY` |
| `champsim_config.json` | Where to set `"prefetcher": "fourier"` for L1D |
