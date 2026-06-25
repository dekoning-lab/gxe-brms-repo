#!/usr/bin/env python3
"""
60_recode_5httlpr_functional.py

Build an alternate analysis-ready dataset in which 5-HTTLPR is coded under
the Hu et al. (2006) triallelic-functional model rather than the biallelic
S-count model.

CODING DIFFERENCE
=================
Original (biallelic S-count) coding in column `5HTTLPR` of
`Fulldata_with_PCs_and_maternal_PCs.txt`:
    0 = no S alleles   (LA/LA, LA/LG, LG/LA, LG/LG)
    1 = one S allele   (S/LA, S/LG)
    2 = two S alleles  (S/S)

New (triallelic-functional) coding under Hu et al. 2006 — both S and LG
are treated as functionally short (rs25531 G in the long allele reduces
transcriptional efficiency to a level comparable to the S allele):
    0 = LA/LA
    1 = S/LA, LG/LA, LA/LG  (one functionally-short allele)
    2 = S/S, S/LG, LG/S, LG/LG  (two functionally-short alleles)

PIPELINE
========
1. Read the Kobor Excel sheet
   (`data/Kobor_lab_genotype_summary_5HTTLPR_triallelic.xlsx`, Sheet1).
   Sample IDs are in column "Samples"; HTTLPR genotypes in column 9.
2. Parse each triallelic genotype string ("S/La", "S/Lg", "La/Lg", etc.)
   into an unordered pair of allele tokens {S, LA, LG}.
3. Compute both (a) the old biallelic S-count (for sanity-check against
   the existing dataset) and (b) the new functional-S count.
4. Read `data/Fulldata_with_PCs_and_maternal_PCs.txt`, join by sample ID.
5. For rows with no LG allele in the triallelic call, the new functional
   coding must equal the old biallelic S-count exactly — verify and abort
   on any mismatch (these would indicate a join error or upstream
   miscoding).
6. Replace the `5HTTLPR` column in-place with the new functional coding.
7. Write to `data_alt_5httlpr/Fulldata_with_PCs_and_maternal_PCs_alt5httlpr.txt`.
8. Also write a verification report.

USAGE
=====
Run from: the repository root:
    python3 scripts/60_recode_5httlpr_functional.py

OUTPUT
======
- data_alt_5httlpr/Fulldata_with_PCs_and_maternal_PCs_alt5httlpr.txt
- data_alt_5httlpr/recoding_report.txt
- data_alt_5httlpr/recoding_join_audit.csv
"""
from __future__ import annotations
import csv
import os
import sys
from collections import Counter
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT  = ROOT / "data_alt_5httlpr"
OUT.mkdir(parents=True, exist_ok=True)

EXCEL_PATH = DATA / "Kobor_lab_genotype_summary_5HTTLPR_triallelic.xlsx"
INPUT_TXT  = DATA / "Fulldata_with_PCs_and_maternal_PCs.txt"
OUTPUT_TXT = OUT / "Fulldata_with_PCs_and_maternal_PCs_alt5httlpr.txt"
REPORT     = OUT / "recoding_report.txt"
AUDIT_CSV  = OUT / "recoding_join_audit.csv"

# ---- Parsing rules -----------------------------------------------------
ALLELE_TOKENS = {"S", "LA", "LG"}

def parse_genotype(s: str | None) -> tuple[str, str] | None:
    """Parse a triallelic genotype string like 'S/La', 'Lg/La' into a sorted
    pair of tokens drawn from {S, LA, LG}. Returns None for blanks/None.
    Aborts on unrecognised token, so any new variant in the source data
    surfaces as a hard error rather than silent miscoding."""
    if s is None:
        return None
    s = str(s).strip()
    if not s:
        return None
    if "/" not in s:
        raise ValueError(f"Bad genotype (no '/'): {s!r}")
    a, b = s.split("/", 1)
    a, b = a.strip().upper(), b.strip().upper()
    if a not in ALLELE_TOKENS or b not in ALLELE_TOKENS:
        raise ValueError(f"Unknown allele in genotype {s!r}: tokens={a},{b}")
    # Stable canonical order (alphabetical) to collapse e.g. 'La/Lg' == 'Lg/La'
    return tuple(sorted([a, b]))

def count_S_old(pair: tuple[str, str] | None) -> int | None:
    """Old biallelic coding: count of S alleles (LG collapsed to L)."""
    if pair is None:
        return None
    return sum(1 for a in pair if a == "S")

def count_functional_S(pair: tuple[str, str] | None) -> int | None:
    """New Hu-et-al-2006 functional coding: count of {S, LG} alleles."""
    if pair is None:
        return None
    return sum(1 for a in pair if a in ("S", "LG"))

def has_LG(pair: tuple[str, str] | None) -> bool:
    if pair is None:
        return False
    return "LG" in pair

# ---- Load Excel --------------------------------------------------------
def load_excel(path: Path) -> dict[str, str | None]:
    """Return {sample_id: raw_genotype_str_or_None}."""
    wb = openpyxl.load_workbook(str(path), data_only=True)
    ws = wb["Sheet1"]
    rows = list(ws.iter_rows(values_only=True))
    header = rows[0]
    if header[0] != "Samples":
        raise ValueError(f"Unexpected header[0]: {header[0]!r}")
    if not (header[8] and "HTTLPR" in str(header[8])):
        raise ValueError(f"Unexpected header[8]: {header[8]!r}")
    out: dict[str, str | None] = {}
    for r in rows[1:]:
        sid = r[0]
        if sid is None or str(sid).strip() == "":
            continue
        sid = str(sid).strip()
        gt = r[8]
        out[sid] = (str(gt).strip() if gt is not None else None)
    return out

# ---- Load TXT ----------------------------------------------------------
def load_txt(path: Path) -> tuple[list[str], list[list[str]]]:
    with open(path, "r", newline="") as f:
        rdr = csv.reader(f, delimiter="\t")
        rows = list(rdr)
    header, body = rows[0], rows[1:]
    return header, body

# ---- Main --------------------------------------------------------------
def main() -> int:
    print(f"Reading Excel: {EXCEL_PATH}")
    raw_geno = load_excel(EXCEL_PATH)
    print(f"  {len(raw_geno)} samples in Excel")

    # Parse and tabulate
    parsed: dict[str, tuple[str, str] | None] = {}
    raw_counts: Counter[str | None] = Counter()
    for sid, gt in raw_geno.items():
        try:
            p = parse_genotype(gt)
        except ValueError as e:
            print(f"  ERROR parsing sample {sid}: {e}", file=sys.stderr)
            return 1
        parsed[sid] = p
        raw_counts[gt if gt is not None else None] += 1

    print("\nTriallelic genotype frequencies (canonicalised):")
    canon_counts: Counter = Counter()
    for sid, p in parsed.items():
        canon_counts[p] += 1
    for k, n in sorted(canon_counts.items(),
                       key=lambda x: (-x[1] if x[0] is not None else -1)):
        label = "/".join(k) if k is not None else "NA"
        old = count_S_old(k)
        new = count_functional_S(k)
        print(f"  {label:8} n={n:>3}   old_S={old}   new_functional={new}")

    # Sanity-check expected totals (from README)
    expected_old: dict[int, int] = {0: 60, 1: 92, 2: 46}
    expected_new: dict[int, int] = {0: 45, 1: 95, 2: 58}
    got_old: Counter[int | None] = Counter()
    got_new: Counter[int | None] = Counter()
    for p in parsed.values():
        got_old[count_S_old(p)] += 1
        got_new[count_functional_S(p)] += 1
    print("\nDerived counts vs README:")
    print(f"  old biallelic 0={got_old[0]} (expected 60), "
          f"1={got_old[1]} (expected 92), 2={got_old[2]} (expected 46), "
          f"NA={got_old[None]}")
    print(f"  new functional 0={got_new[0]} (expected 45), "
          f"1={got_new[1]} (expected 95), 2={got_new[2]} (expected 58), "
          f"NA={got_new[None]}")
    for k, exp in expected_old.items():
        if got_old[k] != exp:
            print(f"  WARNING: old count[{k}]={got_old[k]} != expected {exp}",
                  file=sys.stderr)
    for k, exp in expected_new.items():
        if got_new[k] != exp:
            print(f"  WARNING: new count[{k}]={got_new[k]} != expected {exp}",
                  file=sys.stderr)

    # ---- Load TXT ------------------------------------------------------
    print(f"\nReading TXT: {INPUT_TXT}")
    header, body = load_txt(INPUT_TXT)
    print(f"  {len(body)} sample rows, {len(header)} columns")

    samples_col = header.index("Samples")
    httl_col    = header.index("5HTTLPR")
    print(f"  Samples column index = {samples_col}")
    print(f"  5HTTLPR column index = {httl_col}")

    # ---- Join + validate ----------------------------------------------
    audit_rows: list[dict] = []
    mismatches: list[tuple[str, int | None, int | None, str | None]] = []
    in_txt_not_excel: list[str] = []
    in_excel_not_txt: list[str] = []
    new_col_values: list[str] = []

    txt_ids = set()
    for r in body:
        sid = r[samples_col].strip()
        txt_ids.add(sid)
        old_val_txt_raw = r[httl_col].strip()  # current value in TXT
        if old_val_txt_raw in ("", "NA"):
            old_val_txt: int | None = None
        else:
            try:
                old_val_txt = int(float(old_val_txt_raw))
            except ValueError:
                old_val_txt = None

        if sid not in parsed:
            in_txt_not_excel.append(sid)
            # No triallelic call available; preserve original missingness.
            # Empty string matches the original TXT's NA convention (R reads it
            # as NA for numeric columns).
            new_col_values.append("")
            audit_rows.append({
                "Samples": sid, "in_excel": False,
                "raw_geno": "", "canon_geno": "",
                "has_LG": False,
                "old_in_txt": old_val_txt if old_val_txt is not None else "NA",
                "old_derived": "",
                "new_functional": "NA",
                "lg_distinct_status": "no_excel_call"
            })
            continue

        pair = parsed[sid]
        old_derived = count_S_old(pair)
        new_val_from_excel = count_functional_S(pair)

        # POLICY: preserve the original analytic sample exactly.
        # If TXT was NA, write NA in the new column too (regardless of what
        # Excel says). This isolates the coding-change effect from any
        # sample-change effect (the new coding is applied only where the
        # original analysis had observed data).
        # Match the original TXT format: empty string for NA, "X.0" floats
        # for integer values (R's read.table parses both).
        def _format_int(v):
            return f"{int(v)}.0"
        if old_val_txt is None:
            new_col_value = ""        # empty == NA in the source TXT format
            status = "txt_was_NA_preserved"
            if pair is not None:
                mismatches.append(
                    (sid, None, old_derived,
                     "TXT NA preserved (Excel had call; not used to preserve sample)")
                )
        else:
            if pair is None:
                new_col_value = ""
                status = "TXT_value_Excel_NA"
                mismatches.append(
                    (sid, old_val_txt, None,
                     "TXT non-NA but Excel NA: cannot recode")
                )
            else:
                new_col_value = _format_int(new_val_from_excel)
                if has_LG(pair):
                    status = "LG_present_recoded"
                else:
                    status = "LG_absent_should_match"
                    if old_derived != old_val_txt:
                        mismatches.append(
                            (sid, old_val_txt, old_derived,
                             "TXT old != Excel-derived old (LG absent)")
                        )

        new_col_values.append(new_col_value)

        audit_rows.append({
            "Samples": sid, "in_excel": True,
            "raw_geno": raw_geno.get(sid, "") or "",
            "canon_geno": "/".join(pair) if pair is not None else "",
            "has_LG": has_LG(pair),
            "old_in_txt": old_val_txt if old_val_txt is not None else "NA",
            "old_derived": old_derived if old_derived is not None else "NA",
            "new_functional": new_col_value,
            "lg_distinct_status": status
        })

    in_excel_not_txt = [s for s in parsed if s not in txt_ids]

    # ---- Report ---------------------------------------------------------
    print(f"\nJoin summary:")
    print(f"  TXT sample IDs: {len(txt_ids)}")
    print(f"  Excel sample IDs: {len(parsed)}")
    print(f"  Matched: {len(txt_ids & set(parsed))}")
    print(f"  In TXT, not in Excel: {len(in_txt_not_excel)}")
    print(f"  In Excel, not in TXT: {len(in_excel_not_txt)}")
    if in_txt_not_excel:
        print("    TXT-only:", in_txt_not_excel[:10], "..." if len(in_txt_not_excel)>10 else "")
    if in_excel_not_txt:
        print("    Excel-only:", in_excel_not_txt[:10], "..." if len(in_excel_not_txt)>10 else "")

    print(f"\nValidation: rows where LG-absent triallelic SHOULD match old TXT coding:")
    n_lg_absent = sum(1 for r in audit_rows
                      if r["lg_distinct_status"] == "LG_absent_should_match")
    n_match_clean = sum(1 for r in audit_rows
                        if r["lg_distinct_status"] == "LG_absent_should_match"
                        and r["old_in_txt"] == r["old_derived"])
    print(f"  total {n_lg_absent}, matching {n_match_clean}, "
          f"mismatches {n_lg_absent - n_match_clean}")
    if mismatches:
        print("\nMISMATCH / NOTE log:")
        for sid, txt, derived, note in mismatches[:20]:
            print(f"  {sid}  txt={txt}  derived={derived}  ({note})")
        if len(mismatches) > 20:
            print(f"  ... and {len(mismatches)-20} more")
        # Only abort on the truly bad cases — TXT-old differs from Excel-derived
        # in an LG-absent row, OR TXT had a value but Excel has none and we
        # cannot recode.
        hard = [m for m in mismatches
                if "TXT old != Excel-derived old (LG absent)" in m[3]
                or "TXT non-NA but Excel NA" in m[3]]
        if hard:
            print("\nABORTING: hard mismatches found:", file=sys.stderr)
            for m in hard:
                print(f"  {m}", file=sys.stderr)
            _write_audit(audit_rows)
            return 2

    def _to_int_or_na(v):
        if v in ("NA", "", None):
            return "NA"
        try:
            return int(float(v))
        except (TypeError, ValueError):
            return "NA"

    # ---- Marginal change accounting -----------------------------------
    deltas = Counter()
    n_changed = 0
    for r in audit_rows:
        o = _to_int_or_na(r["old_in_txt"])
        n = _to_int_or_na(r["new_functional"])
        if o == "NA" or n == "NA":
            continue
        if o != n:
            n_changed += 1
            deltas[(o, n)] += 1
    print(f"\nRows whose 5HTTLPR value changes from old to new: {n_changed}")
    for (o, n), c in sorted(deltas.items()):
        print(f"  {o} → {n}: {c} rows")

    # ---- Frequency comparison summary ---------------------------------
    print("\nFrequency comparison (in TXT rows matched to Excel):")
    old_freq: Counter = Counter()
    new_freq: Counter = Counter()
    for r in audit_rows:
        if not r["in_excel"]:
            continue
        old_freq[_to_int_or_na(r["old_in_txt"])] += 1
        new_freq[_to_int_or_na(r["new_functional"])] += 1
    print(f"  Old:  0={old_freq[0]}  1={old_freq[1]}  2={old_freq[2]}  NA={old_freq['NA']}")
    print(f"  New:  0={new_freq[0]}  1={new_freq[1]}  2={new_freq[2]}  NA={new_freq['NA']}")

    # ---- Write output TXT ---------------------------------------------
    out_body = []
    for r, new_val in zip(body, new_col_values):
        new_r = list(r)
        new_r[httl_col] = new_val
        out_body.append(new_r)

    with open(OUTPUT_TXT, "w", newline="") as f:
        # Use minimal quoting to mirror the source TXT (which is plain TSV).
        w = csv.writer(f, delimiter="\t", lineterminator="\n",
                       quoting=csv.QUOTE_MINIMAL)
        w.writerow(header)
        for r in out_body:
            w.writerow(r)
    print(f"\nWrote {OUTPUT_TXT} ({len(out_body)} rows)")

    # ---- Write audit CSV ---------------------------------------------
    _write_audit(audit_rows)

    # ---- Write report -------------------------------------------------
    with open(REPORT, "w") as f:
        f.write("5-HTTLPR recoding report\n")
        f.write("========================\n\n")
        f.write(f"Source Excel:  {EXCEL_PATH}\n")
        f.write(f"Source TXT:    {INPUT_TXT}\n")
        f.write(f"Output TXT:    {OUTPUT_TXT}\n\n")
        f.write("Coding switch: biallelic S-count -> Hu et al. 2006 triallelic-functional\n")
        f.write("  S and LG both counted as functionally short.\n\n")
        f.write(f"Excel samples:     {len(parsed)}\n")
        f.write(f"TXT samples:       {len(txt_ids)}\n")
        f.write(f"Matched samples:   {len(txt_ids & set(parsed))}\n")
        f.write(f"TXT-only samples:  {len(in_txt_not_excel)}\n")
        f.write(f"Excel-only samples:{len(in_excel_not_txt)}\n\n")
        f.write(f"Rows changed by recoding (LG-bearing genotypes): {n_changed}\n")
        for (o, n), c in sorted(deltas.items()):
            f.write(f"  {o} -> {n}: {c} rows\n")
        f.write("\nFrequency comparison in matched rows:\n")
        f.write(f"  old:  0={old_freq[0]}  1={old_freq[1]}  2={old_freq[2]}  NA={old_freq['NA']}\n")
        f.write(f"  new:  0={new_freq[0]}  1={new_freq[1]}  2={new_freq[2]}  NA={new_freq['NA']}\n")
    print(f"Wrote {REPORT}")

    return 0

def _write_audit(audit_rows):
    with open(AUDIT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(audit_rows[0].keys()))
        w.writeheader()
        for r in audit_rows:
            w.writerow(r)
    print(f"Wrote {AUDIT_CSV}")

if __name__ == "__main__":
    sys.exit(main())
