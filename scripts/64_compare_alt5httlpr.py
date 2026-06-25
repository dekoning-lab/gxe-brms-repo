#!/usr/bin/env python3
"""
64_compare_alt5httlpr.py

Compare draw-pooled posterior summaries under the production biallelic
S-count coding of 5-HTTLPR vs the alternate Hu et al. (2006)
triallelic-functional coding. Reports deltas across all parameters and
focuses on the PP >= 70% headline effects.

Inputs:
  results/v3_m4/<pv>/summary.csv               (original)
  results_alt_5httlpr/v3_m4/<pv>/summary.csv   (alternate)

Output:
  results_alt_5httlpr/comparison_original_vs_alt5httlpr.csv
  results_alt_5httlpr/comparison_report.txt
"""
import csv
import os
import sys

ORIG_DIR = "results/v3_m4"
ALT_DIR  = "results_alt_5httlpr/v3_m4"
OUT_CSV  = "results_alt_5httlpr/comparison_original_vs_alt5httlpr.csv"
PVS = ["sens", "cont", "unre"]
PP_THRESHOLD = 70.0

def load(path):
    rows = {}
    with open(path) as f:
        rdr = csv.DictReader(f)
        for r in rdr:
            rows[r["term"]] = r
    return rows

all_rows = []
for pv in PVS:
    orig = load(os.path.join(ORIG_DIR, pv, "summary.csv"))
    alt_path = os.path.join(ALT_DIR, pv, "summary.csv")
    if not os.path.exists(alt_path):
        print(f"  MISSING: {alt_path}", file=sys.stderr)
        continue
    alt = load(alt_path)
    terms = list(orig.keys())
    for t in terms:
        if t not in alt:
            continue
        o = orig[t]; a = alt[t]
        try:
            o_est = float(o["estimate"]); a_est = float(a["estimate"])
        except ValueError:
            continue
        d_est = a_est - o_est
        try:
            o_pp = 100*float(o["post_prob"]) if o["post_prob"] not in ("","NA") else None
            a_pp = 100*float(a["post_prob"]) if a["post_prob"] not in ("","NA") else None
        except ValueError:
            o_pp = a_pp = None
        d_pp = (a_pp - o_pp) if (o_pp is not None and a_pp is not None) else None
        flip = (o.get("direction","") != a.get("direction","")) if o.get("direction") and a.get("direction") else False
        try:
            o_q025 = float(o["q025"]); o_q975 = float(o["q975"])
            a_q025 = float(a["q025"]); a_q975 = float(a["q975"])
        except ValueError:
            o_q025 = o_q975 = a_q025 = a_q975 = None
        all_rows.append({
            "parenting": pv,
            "term": t,
            "effect_type": o.get("effect_type",""),
            "orig_estimate": o_est, "alt_estimate": a_est,
            "delta_estimate": d_est,
            "orig_q025": o_q025, "orig_q975": o_q975,
            "alt_q025": a_q025, "alt_q975": a_q975,
            "orig_post_prob_pct": o_pp, "alt_post_prob_pct": a_pp,
            "delta_post_prob_pp": d_pp,
            "orig_direction": o.get("direction",""), "alt_direction": a.get("direction",""),
            "direction_flip": flip,
        })

os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
with open(OUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
    w.writeheader()
    for r in all_rows: w.writerow(r)
print(f"Wrote {OUT_CSV}: {len(all_rows)} rows\n")

print("="*80)
print("Top 10 by |Δ posterior mean| (across all parameters)")
print("="*80)
all_e = [r for r in all_rows if r.get("delta_estimate") is not None]
for r in sorted(all_e, key=lambda x: abs(x["delta_estimate"]), reverse=True)[:10]:
    print(f"  {r['parenting']:5} | {r['term']:24} | orig={r['orig_estimate']:+.4f} "
          f"alt={r['alt_estimate']:+.4f} Δ={r['delta_estimate']:+.4f} "
          f"PP {r['orig_post_prob_pct']:5.1f}→{r['alt_post_prob_pct']:5.1f}")

print()
print("="*80)
print("Top 10 by |Δ posterior probability| (percentage points)")
print("="*80)
all_p = [r for r in all_rows if r.get("delta_post_prob_pp") is not None]
for r in sorted(all_p, key=lambda x: abs(x["delta_post_prob_pp"]), reverse=True)[:10]:
    print(f"  {r['parenting']:5} | {r['term']:24} | orig={r['orig_estimate']:+.4f} "
          f"alt={r['alt_estimate']:+.4f} Δ={r['delta_estimate']:+.4f} "
          f"PP {r['orig_post_prob_pct']:5.1f}→{r['alt_post_prob_pct']:5.1f} ΔPP={r['delta_post_prob_pp']:+.2f}")

print()
print("="*80)
print(f"Effects with PP >= {PP_THRESHOLD}% in EITHER pool — change under alt coding")
print("="*80)
header = f"  {'pv':5} | {'term':24} | {'orig est (PP%)':14} | {'alt est (PP%)':14} | {'Δest':>7} | {'ΔPP':>5}"
print(header)
print("  " + "-"*(len(header)-2))
for r in all_rows:
    o_pp = r["orig_post_prob_pct"]; a_pp = r["alt_post_prob_pct"]
    if (o_pp is not None and o_pp >= PP_THRESHOLD) or (a_pp is not None and a_pp >= PP_THRESHOLD):
        flag = " *FLIP*" if r["direction_flip"] else ""
        # mark if effect involves X5HTTLPR (the locus whose coding changed)
        if "X5HTTLPR" in r["term"]:
            flag += "  [5-HTTLPR effect]"
        print(f"  {r['parenting']:5} | {r['term']:24} | "
              f"{r['orig_estimate']:+.3f} ({o_pp:4.1f})  | "
              f"{r['alt_estimate']:+.3f} ({a_pp:4.1f})  | "
              f"{r['delta_estimate']:+7.4f} | {r['delta_post_prob_pp']:+5.2f}{flag}")

print()
print("="*80)
print("All 5-HTTLPR effects (direct effect of the coding change)")
print("="*80)
print(f"  {'pv':5} | {'term':24} | {'orig est (PP%)':14} | {'alt est (PP%)':14} | {'Δest':>7} | {'ΔPP':>5}")
print("  " + "-"*70)
for r in all_rows:
    if "X5HTTLPR" not in r["term"]:
        continue
    print(f"  {r['parenting']:5} | {r['term']:24} | "
          f"{r['orig_estimate']:+.3f} ({r['orig_post_prob_pct']:4.1f})  | "
          f"{r['alt_estimate']:+.3f} ({r['alt_post_prob_pct']:4.1f})  | "
          f"{r['delta_estimate']:+7.4f} | {r['delta_post_prob_pp']:+5.2f}")

flips = [r for r in all_rows if r["direction_flip"]]
print(f"\nDirection flips: {len(flips)}")
for r in flips:
    print(f"  {r['parenting']:5} | {r['term']:24} | "
          f"orig {r['orig_direction']} ({r['orig_post_prob_pct']:5.2f}%)  →  "
          f"alt  {r['alt_direction']} ({r['alt_post_prob_pct']:5.2f}%)")
