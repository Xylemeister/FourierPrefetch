# The Goertzel Algorithm

## The problem it solves

A full DFT takes a sequence of N samples and computes power at **all** N frequency bins. That's O(N^2) work (or O(N log N) with FFT). But we don't care about all frequencies -- we only want to check a few specific periods (2, 3, 4, 6, 8, 12). Goertzel computes **one frequency bin** in O(N) with just 1 multiply-add per sample.

## How it works

Say you have a buffer of N=24 deltas and you want to check if there's a repeating pattern with period P=4. That corresponds to frequency k = N/P = 6.

The algorithm uses a single precomputed coefficient:
```
coeff = 2 * cos(2 * pi * k / N)
```

Then it runs a simple recurrence over all N samples:
```
s0 = 0, s1 = 0

for each sample x[i]:
    s_new = x[i] + coeff * s0 - s1
    s1 = s0
    s0 = s_new

power = s0^2 + s1^2 - coeff * s0 * s1
```

That's it -- just N multiply-adds and a final power computation. The output `power` tells you how much energy the signal has at that frequency.

## What the power value means

- **High power at period P** means the deltas repeat every P accesses (e.g., `[4, 8, 4, 8]` has high power at P=2)
- **Low power everywhere** means no periodic pattern, deltas are random or non-repeating
- **High power at multiple periods** means the highest one is the fundamental period

## Why it works for prefetching

Consider a nested loop that produces this delta sequence:
```
[1, 2, 3, 1, 2, 3, 1, 2, 3]
```

A simple stride detector sees deltas changing every access and gives up. Goertzel at period P=3 would return high power, telling us "this signal repeats every 3 steps." We then extract the last 3 deltas `[1, 2, 3]` and replay them cyclically to predict the next access.

## Our specific setup

We test 6 candidate periods with precomputed coefficients:

| Period P | Frequency k = 24/P | coeff = 2*cos(2*pi*k/24) |
|----------|--------------------:|--------------------------|
| 2        | 12                  | -2.0                     |
| 3        | 8                   | -1.0                     |
| 4        | 6                   | 0.0                      |
| 6        | 4                   | 1.0                      |
| 8        | 3                   | 1.4142                   |
| 12       | 2                   | 1.7321                   |

Total cost: 6 periods x 24 samples = 144 multiply-adds per analysis. Negligible compared to the cost of simulating a single cache access in ChampSim.

## C++ implementation

```cpp
// For buffer of N samples, test frequency k = N/P
// coeff = 2 * cos(2 * pi * k / N), precomputed per period
float goertzel_power(const int16_t* buf, int head, int N, float coeff) {
    float s0 = 0, s1 = 0;
    for (int i = 0; i < N; i++) {
        float s_new = buf[(head + i) % N] + coeff * s0 - s1;
        s1 = s0;
        s0 = s_new;
    }
    return s0*s0 + s1*s1 - coeff*s0*s1;  // |X(k)|^2
}
```
