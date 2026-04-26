#!/usr/bin/env python3
"""
Parse [ISSUED OFFSET] lines from a ChampSim log and plot offset frequency.

Usage:
    ./bin/... 2>&1 | python3 plot_offsets.py
    python3 plot_offsets.py output.txt
    python3 plot_offsets.py results/          # aggregates all .txt files recursively
    python3 plot_offsets.py output.txt --bins 100 --out offsets.png
"""

import sys
import re
import argparse
import collections
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # no display required
import matplotlib.pyplot as plt

FIGURE_DIR = Path(__file__).parent / 'figure' / 'offsets'


def parse_offsets(lines):
    pattern = re.compile(r'\[ISSUED OFFSET\]\s+(-?\d+)')
    offsets = []
    for line in lines:
        m = pattern.search(line)
        if m:
            offsets.append(int(m.group(1)))
    return offsets


def collect_offsets(path):
    p = Path(path)
    if p.is_file():
        files = [p]
    elif p.is_dir():
        files = sorted(p.rglob('*.txt'))
        if not files:
            print(f"No .txt files found under {p}", file=sys.stderr)
            sys.exit(1)
        print(f"Found {len(files)} file(s) under {p}")
    else:
        print(f"'{path}' is not a file or directory", file=sys.stderr)
        sys.exit(1)

    offsets = []
    for f in files:
        with open(f) as fh:
            chunk = parse_offsets(fh)
        print(f"  {f.name}: {len(chunk)} offset(s)")
        offsets.extend(chunk)
    return offsets


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('path', nargs='?', help='Log file or directory (default: stdin)')
    parser.add_argument('--bins', type=int, default=None,
                        help='Number of histogram bins (default: auto)')
    parser.add_argument('--out', default=None,
                        help='Override output filename (still saved inside figure/offsets/)')
    parser.add_argument('--top', type=int, default=None,
                        help='Only plot the N most frequent offset values')
    args = parser.parse_args()

    if args.path:
        offsets = collect_offsets(args.path)
    else:
        offsets = parse_offsets(sys.stdin)
    if not offsets:
        print("No [ISSUED OFFSET] lines found.", file=sys.stderr)
        sys.exit(1)

    print(f"Total offsets parsed: {len(offsets)}")
    print(f"Range: [{min(offsets)}, {max(offsets)}]")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # --- Left: histogram over the raw offset values ---
    ax = axes[0]
    n_bins = args.bins or min(200, len(set(offsets)))
    ax.hist(offsets, bins=n_bins, color='steelblue', edgecolor='none')
    ax.set_xlabel('Issued offset (cache lines)')
    ax.set_ylabel('Frequency')
    ax.set_title('Offset distribution (histogram)')
    ax.grid(axis='y', linewidth=0.4, alpha=0.6)

    # --- Right: top-N most frequent discrete values ---
    ax2 = axes[1]
    counter = collections.Counter(offsets)
    top_n = args.top or min(40, len(counter))
    most_common = counter.most_common(top_n)
    most_common.sort(key=lambda x: x[0])          # sort by offset value
    vals, freqs = zip(*most_common)
    ax2.bar(range(len(vals)), freqs, color='darkorange', edgecolor='none')
    ax2.set_xticks(range(len(vals)))
    ax2.set_xticklabels([str(v) for v in vals], rotation=90, fontsize=7)
    ax2.set_xlabel('Offset value')
    ax2.set_ylabel('Count')
    ax2.set_title(f'Top {top_n} most frequent offsets')
    ax2.grid(axis='y', linewidth=0.4, alpha=0.6)

    fig.tight_layout()

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    if args.out:
        out_path = FIGURE_DIR / args.out
    else:
        stem = Path(args.path).stem if args.path else 'stdin'
        out_path = FIGURE_DIR / f'{stem}_offsets.png'

    fig.savefig(out_path, dpi=150)
    print(f"Saved to {out_path}")


if __name__ == '__main__':
    main()
