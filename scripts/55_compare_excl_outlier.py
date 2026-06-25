#!/usr/bin/env python3
"""
Compare original draw-pooled estimates (using all 100 imputations) to the
outlier-excluded estimates (99 imputations; sens drops imp 3, cont drops
imp 1, unre drops imp 38).

For each parameter:
  - delta_est       = excluded_estimate - original_estimate (log-odds)
  - delta_post_prob = (excluded_post_prob - original_post_prob) * 100  (pp)
  - direction_flip  = original.direction != excluded.direction (T/F)

Writes one combined CSV and prints a compact summary of the largest
deltas and any direction flips.

Run from: the repository root
Usage: python3 scripts/55_compare_excl_outlier.py
"""
import csv
import os
import sys

ORIG_DIR = "results/v3_m4"
EXCL_DIR = "results/v3_m4_excl_outlier"
OUT_CSV  = "results/v3_m4_excl_outlier/comparison_original_vs_excluded.csv"

PVS = ["sens", "cont", "unre"]
EXCLUDED = {"sens": 3, "cont": 1, "unre": 38}
PP_THRESHOLD = 70.0  # percent

def load(path):
    rows = {}
    with open(path) as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            rows[r["term"]] = r
    return rows

def fnum(x, nd=4):
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return ""

all_rows = []
for pv in PVS:
    orig_path = os.path.join(ORIG_DIR, pv, "summary.csv")
    excl_path = os.path.join(EXCL_DIR, pv, "summary.csv")
    if not os.path.exists(excl_path):
        print(f"  MISSING: {excl_path}", file=sys.stderr)
        continue
    orig = load(orig_path)
    excl = load(excl_path)
    terms = list(orig.keys())
    for t in terms:
        if t not in excl:
            continue
        o = orig[t]
        e = excl[t]
        try:
            o_est = float(o["estimate"]); e_est = float(e["estimate"])
        except ValueError:
            continue
        d_est = e_est - o_est
        try:
            o_pp = 100*float(o["post_prob"]) if o["post_prob"] not in ("", "NA") else None
            e_pp = 100*float(e["post_prob"]) if e["post_prob"] not in ("", "NA") else None
        except ValueError:
            o_pp = e_pp = None
        d_pp = (e_pp - o_pp) if (o_pp is not None and e_pp is not None) else None
        flip = (o.get("direction", "") != e.get("direction", "")) if o.get("direction") and e.get("direction") else False
        try:
            o_q025 = float(o["q025"]); o_q975 = float(o["q975"])
            e_q025 = float(e["q025"]); e_q975 = float(e["q975"])
        except ValueError:
            o_q025 = o_q975 = e_q025 = e_q975 = None
        all_rows.append({
            "parenting": pv,
            "excluded_imp": EXCLUDED[pv],
            "term": t,
            "effect_type": o.get("effect_type", ""),
            "orig_estimate": o_est,
            "excl_estimate": e_est,
            "delta_estimate": d_est,
            "orig_q025": o_q025, "orig_q975": o_q975,
            "excl_q025": e_q025, "excl_q975": e_q975,
            "orig_post_prob_pct": o_pp,
            "excl_post_prob_pct": e_pp,
            "delta_post_prob_pp": d_pp,
            "orig_direction": o.get("direction", ""),
            "excl_direction": e.get("direction", ""),
            "direction_flip": flip,
        })

os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
with open(OUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
    w.writeheader()
    for r in all_rows:
        w.writerow(r)
print(f"Wrote {OUT_CSV}: {len(all_rows)} rows")

# ── Compact summary ─────────────────────────────────────────────
print()
print("=" * 80)
print("Max absolute deltas across all parameters")
print("=" * 80)
all_d_est = [r for r in all_rows if r.get("delta_estimate") is not None]
all_d_pp = [r for r in all_rows if r.get("delta_post_prob_pp") is not None]
top_est = sorted(all_d_est, key=lambda r: abs(r["delta_estimate"]), reverse=True)[:8]
top_pp  = sorted(all_d_pp,  key=lambda r: abs(r["delta_post_prob_pp"]), reverse=True)[:8]
print("\nTop 8 by |delta_estimate| (log-odds):")
for r in top_est:
    print(f"  {r['parenting']:5} | {r['term']:24} | orig={r['orig_estimate']:+.4f}  "
          f"excl={r['excl_estimate']:+.4f}  delta={r['delta_estimate']:+.4f}  "
          f"PP {r['orig_post_prob_pct']:5.1f}→{r['excl_post_prob_pct']:5.1f}")
print("\nTop 8 by |delta PP| (percentage points):")
for r in top_pp:
    print(f"  {r['parenting']:5} | {r['term']:24} | orig={r['orig_estimate']:+.4f}  "
          f"excl={r['excl_estimate']:+.4f}  delta={r['delta_estimate']:+.4f}  "
          f"PP {r['orig_post_prob_pct']:5.1f}→{r['excl_post_prob_pct']:5.1f}  "
          f"dPP={r['delta_post_prob_pp']:+.2f}pp")

# ── Effects with PP >= 70 at the primary scale ─────────────────
print("\n" + "=" * 80)
print(f"Effects with PP >= {PP_THRESHOLD}% in the ORIGINAL pool — change under exclusion")
print("=" * 80)
header = f"  {'pv':5} | {'term':24} | {'orig est (PP%)':14} | {'excl est (PP%)':14} | {'Δest':>7} | {'ΔPP':>5}"
print(header)
print("  " + "-" * (len(header) - 2))
for r in all_rows:
    if r["orig_post_prob_pct"] is not None and r["orig_post_prob_pct"] >= PP_THRESHOLD:
        flip_marker = " *FLIP*" if r["direction_flip"] else ""
        print(f"  {r['parenting']:5} | {r['term']:24} | "
              f"{r['orig_estimate']:+.3f} ({r['orig_post_prob_pct']:4.1f})  | "
              f"{r['excl_estimate']:+.3f} ({r['excl_post_prob_pct']:4.1f})  | "
              f"{r['delta_estimate']:+7.4f} | {r['delta_post_prob_pp']:+5.2f}{flip_marker}")

# Direction flips
flips = [r for r in all_rows if r["direction_flip"]]
print()
print(f"Direction flips (any PP): {len(flips)}")
for r in flips:
    print(f"  {r['parenting']:5} | {r['term']:24} | "
          f"orig {r['orig_direction']} ({r['orig_post_prob_pct']:5.2f}%)  →  "
          f"excl {r['excl_direction']} ({r['excl_post_prob_pct']:5.2f}%)")
