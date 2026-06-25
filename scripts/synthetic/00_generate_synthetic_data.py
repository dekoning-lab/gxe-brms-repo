#!/usr/bin/env python3
"""
00_generate_synthetic_data.py

Generate a synthetic dataset matching the schema of the real APrON
analysis file, for the purpose of demonstrating the pipeline without
revealing any real participant data.

DESIGN
------
For each variable in the real input file:
  - Categorical / integer-coded variables (Sex, ethnicity, genetic
    variants, outcome): sample independently from the empirical
    marginal distribution (with proportional missing-rate).
  - Continuous variables (parenting scales, infant age, ancestry PCs):
    sample independently from a smoothed empirical distribution (a
    Gaussian kernel-based draw using the empirical mean + SD).
  - Participant IDs: synthetic `DUMMY-NNNNN` (no overlap with real
    cohort IDs).
  - Maternal IDs: synthetic `MDUMMY-NNNNN`.

The synthetic dataset preserves the input file's column order,
row count (n=210 by default), marginal distributions, and overall
missingness pattern (variable-by-variable). It does NOT preserve
joint dependence between variables — any pairwise correlations
present in the real data are intentionally destroyed.

This means the synthetic data is suitable for:
  - Running the pipeline end-to-end as a smoke test
  - Verifying that file formats, column names, and types parse correctly
  - Adapting the code to other cohorts

It is NOT suitable for:
  - Drawing scientific conclusions
  - Estimating power
  - Anything that depends on the real joint distribution

USAGE
-----
Run from the repository root:

    python3 scripts/synthetic/00_generate_synthetic_data.py

Optional environment variables:
  SYNTHETIC_N          - number of synthetic rows (default 210)
  SYNTHETIC_SEED       - RNG seed (default 20260511)
  SYNTHETIC_REAL_PATH  - path to the real input file (only used as a
                         marginal-distribution template). If absent,
                         the script falls back to a hard-coded
                         "schema-only" template that produces a
                         dataset with the correct columns and types
                         but draws from naive marginal distributions
                         documented inline.

OUTPUT
------
  data/synthetic/Fulldata_synthetic.txt
  data/synthetic/synthetic_data_provenance.txt
"""
from __future__ import annotations
import csv
import math
import os
import sys
import random
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "data" / "synthetic"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_TXT = OUT_DIR / "Fulldata_synthetic.txt"
PROV    = OUT_DIR / "synthetic_data_provenance.txt"

REAL_PATH = Path(
    os.environ.get(
        "SYNTHETIC_REAL_PATH",
        str(ROOT / "data" / "Fulldata_with_PCs_and_maternal_PCs.txt")
    )
)
N_SYNTH = int(os.environ.get("SYNTHETIC_N", "210"))
SEED    = int(os.environ.get("SYNTHETIC_SEED", "20260511"))

random.seed(SEED)

# ---- Column-type registry ---------------------------------------------
# These are the columns the pipeline actually uses. Anything else in the
# input file is reproduced with schema-appropriate dummy values.

CATEGORICAL_INT_COLS = {
    "Sex": {0, 1},
    "ethnicity": None,            # any integer category present in real data
    "Final7PointLikertScale": set(range(1, 8)),
    "AgesABCD": None,             # integer ages
    "Method": None,               # categorical method code
    "FINAL_CLASSIFICATION_(ABCDModel)": None,
    "Score": None,
}

GENETIC_COLS = [
    "DRD2", "BDNF", "CNR1-77", "CNR1-10", "DRD4", "MAOA",
    "SLC6A3-9R", "SLC6A3-10R", "5HTTLPR", "CNR-1-10-A-dom", "CNR-1-10-A-rec",
]

CONTINUOUS_COLS = [
    "infant_age",
    # CBCL scales (raw and T-scores)
    "Emotion_Reactiveraw", "Emotion_ReactiveTscore",
    "Anx_depressedraw", "Anx_depressedTscore",
    "Somatic_Complaintsraw", "Somatic_ComplaintsTscore",
    "Withdrawnraw", "WithdrawnTscore",
    "Sleep_probraw", "Sleep_probTscore",
    "Attn_probraw", "Attn_probTscore",
    "Aggressive_behraw", "Aggressive_behTscore",
    "Internal_probraw", "Internal_probTscore",
    "External_probraw", "External_probTscore",
    "Total_probraw", "Total_probTscore",
    # CARE-Index sub-scales
    "SensCues", "ResDist", "SocFost", "CognFost", "ClarCues", "ResCare",
    "CareTot", "ChildTot", "NCASTtot",
    # CARE-Index primary scales (the analysis variables)
    "sync", "sens", "cont", "unre", "coop", "ccc", "diff", "pass",
]

# All PCs (child + maternal) are continuous
def is_pc_col(c): return c.startswith("PC") or c.startswith("maternal_PC")

ID_COLS = {"Samples", "Inf_APrONID", "MaternalID"}


def empirical_stats(real_path: Path):
    """Read the real file (if available) and return marginal stats per column."""
    if not real_path.exists():
        return None
    with open(real_path) as f:
        rdr = csv.reader(f, delimiter="\t")
        header = next(rdr)
        cols = {h: [] for h in header}
        for row in rdr:
            for h, v in zip(header, row):
                cols[h].append(v)
    return header, cols


def parse_float(s: str):
    if s is None: return None
    s = s.strip()
    if s == "" or s == "NA":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def parse_int(s: str):
    f = parse_float(s)
    if f is None: return None
    return int(f)


def sample_from_floats(values, n, seed_offset):
    """Sample n values from a continuous variable's marginal distribution.

    Approach: rounded-bootstrap with Gaussian smoothing using the empirical
    mean and SD. Preserves the marginal mean and SD exactly in expectation;
    does not preserve higher moments or the exact shape. Falls back to
    constant 0.0 if the column has no observed values."""
    rng = random.Random(SEED + seed_offset)
    miss_p = sum(1 for v in values if parse_float(v) is None) / max(len(values), 1)
    obs = [parse_float(v) for v in values if parse_float(v) is not None]
    if len(obs) == 0:
        return ["" for _ in range(n)]
    mu = sum(obs) / len(obs)
    var = sum((x - mu) ** 2 for x in obs) / max(len(obs) - 1, 1)
    sd = math.sqrt(var) if var > 0 else 0.0
    out = []
    for _ in range(n):
        if rng.random() < miss_p:
            out.append("")
        else:
            v = rng.gauss(mu, sd) if sd > 0 else mu
            out.append(f"{v:.6f}")
    return out


def sample_from_ints(values, n, seed_offset):
    rng = random.Random(SEED + seed_offset)
    miss_p = sum(1 for v in values if parse_int(v) is None) / max(len(values), 1)
    obs = [parse_int(v) for v in values if parse_int(v) is not None]
    if len(obs) == 0:
        return ["" for _ in range(n)]
    out = []
    for _ in range(n):
        if rng.random() < miss_p:
            out.append("")
        else:
            x = rng.choice(obs)  # bootstrap from empirical
            out.append(f"{x}.0")
    return out


# ---- Fallback synthesizer if no real file present ---------------------

def fallback_template(n):
    """Hard-coded schema with documented dummy distributions, used when
    the real file is not available on disk (which is the case for anyone
    cloning the repo without the private data). Produces a dataset with
    the right columns, types, and approximate marginal scales but no
    connection to any real participants."""

    columns = (
        ["Samples", "Sex", "Inf_APrONID", "MaternalID", "infant_age"]
        + ["Emotion_Reactiveraw","Emotion_ReactiveTscore",
           "Anx_depressedraw","Anx_depressedTscore",
           "Somatic_Complaintsraw","Somatic_ComplaintsTscore",
           "Withdrawnraw","WithdrawnTscore",
           "Sleep_probraw","Sleep_probTscore",
           "Attn_probraw","Attn_probTscore",
           "Aggressive_behraw","Aggressive_behTscore",
           "Internal_probraw","Internal_probTscore",
           "External_probraw","External_probTscore",
           "Total_probraw","Total_probTscore",
           "ethnicity"]
        + ["SensCues","ResDist","SocFost","CognFost","ClarCues","ResCare",
           "CareTot","ChildTot","NCASTtot",
           "sync","sens","cont","unre","coop","ccc","diff","pass","Score"]
        + ["DRD2","BDNF","CNR1-77","CNR1-10","DRD4","MAOA",
           "SLC6A3-9R","SLC6A3-10R","5HTTLPR","CNR-1-10-A-dom","CNR-1-10-A-rec"]
        + ["AgesABCD","FINAL_CLASSIFICATION_(ABCDModel)",
           "Final7PointLikertScale","Method"]
        + [f"PC{i}" for i in range(1, 51)]
        + [f"maternal_PC{i}" for i in range(1, 4)]
    )

    rng = random.Random(SEED)
    rows = []
    for i in range(1, n + 1):
        row = {c: "" for c in columns}
        row["Samples"]      = f"DUMMY-{i:05d}"
        row["Inf_APrONID"]  = f"DUMMY-INF-{i:05d}"
        row["MaternalID"]   = f"MDUMMY-{i:05d}"
        row["Sex"]          = f"{rng.choice([0, 1])}.0"
        row["ethnicity"]    = f"{rng.randint(0, 3)}.0"
        row["infant_age"]   = f"{rng.gauss(36.0, 1.5):.4f}"

        # CBCL scales: raw 0–30ish, T-scores around 50 (SD 10)
        for c in ["Emotion_Reactiveraw","Anx_depressedraw","Somatic_Complaintsraw",
                  "Withdrawnraw","Sleep_probraw","Attn_probraw","Aggressive_behraw"]:
            row[c] = f"{max(0, rng.gauss(4, 4)):.2f}"
        for c in ["Emotion_ReactiveTscore","Anx_depressedTscore","Somatic_ComplaintsTscore",
                  "WithdrawnTscore","Sleep_probTscore","Attn_probTscore","Aggressive_behTscore",
                  "Internal_probTscore","External_probTscore","Total_probTscore"]:
            row[c] = f"{rng.gauss(50, 10):.2f}"
        for c in ["Internal_probraw","External_probraw","Total_probraw"]:
            row[c] = f"{max(0, rng.gauss(15, 8)):.2f}"

        # CARE-Index sub-scales: integer 0–14 each
        for c in ["SensCues","ResDist","SocFost","CognFost","ClarCues","ResCare",
                  "CareTot","ChildTot","NCASTtot"]:
            row[c] = f"{rng.randint(0, 14)}.0"

        # CARE-Index primary scales: integer 0–14 each
        for c in ["sync","sens","cont","unre","coop","ccc","diff","pass"]:
            row[c] = f"{rng.randint(0, 14)}.0"
        row["Score"] = f"{rng.gauss(0, 1):.4f}"

        # Genetic variants: 0/1 for most; 0/1/2 for 5HTTLPR
        for c in ["DRD2","BDNF","CNR1-77","CNR1-10","DRD4","MAOA",
                  "SLC6A3-9R","SLC6A3-10R","CNR-1-10-A-dom","CNR-1-10-A-rec"]:
            row[c] = f"{rng.choice([0,0,0,1,1])}.0"
        row["5HTTLPR"] = f"{rng.choice([0,0,1,1,1,2,2])}.0"

        # Outcome and attachment classifications
        row["Final7PointLikertScale"] = f"{rng.randint(1, 7)}.0"
        row["FINAL_CLASSIFICATION_(ABCDModel)"] = f"{rng.randint(1, 4)}.0"
        row["AgesABCD"] = f"{rng.gauss(36, 1.5):.2f}"
        row["Method"]   = "1"

        # Ancestry PCs: standard normal
        for j in range(1, 51):
            row[f"PC{j}"] = f"{rng.gauss(0, 1):.6f}"
        for j in range(1, 4):
            row[f"maternal_PC{j}"] = f"{rng.gauss(0, 1):.6f}"

        # Inject realistic missingness rates (per docs/DATA.md)
        miss_rates = {
            "Final7PointLikertScale": 0.18, "infant_age": 0.15,
            "ethnicity": 0.005,
            "sens": 0.005, "cont": 0.005, "unre": 0.005,
            "DRD2": 0.02, "BDNF": 0.02, "CNR1-77": 0.02,
            "CNR1-10": 0.02, "DRD4": 0.06, "MAOA": 0.03,
            "SLC6A3-9R": 0.025, "SLC6A3-10R": 0.025,
            "5HTTLPR": 0.06,
        }
        for c in [f"PC{i}" for i in range(1,4)]: miss_rates[c] = 0.07
        for c in [f"maternal_PC{i}" for i in range(1,4)]: miss_rates[c] = 0.06
        for c, p in miss_rates.items():
            if rng.random() < p:
                row[c] = ""

        rows.append(row)
    return columns, rows


# ---- Main -------------------------------------------------------------

def main() -> int:
    print(f"Synthetic-data generator")
    print(f"  output: {OUT_TXT}")
    print(f"  n = {N_SYNTH}, seed = {SEED}")
    print(f"  real file template path: {REAL_PATH} "
          f"({'present' if REAL_PATH.exists() else 'ABSENT — using fallback schema'})")

    if REAL_PATH.exists():
        header, cols = empirical_stats(REAL_PATH)
        n_real = len(next(iter(cols.values())))
        print(f"  read {n_real} real rows across {len(header)} columns; using marginals only")

        # Build per-column synthetic vectors
        synth = {h: [] for h in header}
        for idx, h in enumerate(header):
            values = cols[h]
            if h in ID_COLS:
                synth[h] = [
                    f"DUMMY-{i:05d}" if h == "Samples"
                    else f"DUMMY-INF-{i:05d}" if h == "Inf_APrONID"
                    else f"MDUMMY-{i:05d}"
                    for i in range(1, N_SYNTH + 1)
                ]
            elif h in CATEGORICAL_INT_COLS or h in GENETIC_COLS:
                synth[h] = sample_from_ints(values, N_SYNTH, idx)
            elif h in CONTINUOUS_COLS or is_pc_col(h):
                synth[h] = sample_from_floats(values, N_SYNTH, idx)
            else:
                # Unknown column type: treat as float if anything looks numeric,
                # else fall back to constant blank.
                if any(parse_float(v) is not None for v in values):
                    synth[h] = sample_from_floats(values, N_SYNTH, idx)
                else:
                    synth[h] = ["" for _ in range(N_SYNTH)]

        out_rows = []
        for i in range(N_SYNTH):
            out_rows.append([synth[h][i] for h in header])
        columns = header
    else:
        print("  using fallback schema-only template (no marginals from real data)")
        columns, dicts = fallback_template(N_SYNTH)
        out_rows = [[d.get(c, "") for c in columns] for d in dicts]

    # Write output TSV
    with open(OUT_TXT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n",
                       quoting=csv.QUOTE_MINIMAL)
        w.writerow(columns)
        for r in out_rows:
            w.writerow(r)
    print(f"  wrote {OUT_TXT} ({len(out_rows)} rows × {len(columns)} columns)")

    # Provenance record
    with open(PROV, "w") as f:
        f.write("Synthetic dataset provenance\n")
        f.write("============================\n\n")
        f.write(f"Generated by: scripts/synthetic/00_generate_synthetic_data.py\n")
        f.write(f"Seed:         {SEED}\n")
        f.write(f"N rows:       {N_SYNTH}\n")
        real_provided = ("SYNTHETIC_REAL_PATH" in os.environ) or REAL_PATH.exists()
        if real_provided:
            f.write("Real-data template:  (provided via SYNTHETIC_REAL_PATH; "
                    "absolute path not recorded)\n")
        else:
            f.write("Real-data template:  (none; schema-only fallback template used)\n")
        f.write(f"Template was present: {REAL_PATH.exists()}\n\n")
        f.write("This dataset uses synthetic DUMMY-NNNNN participant IDs and\n")
        f.write("contains no real participant data. Marginal distributions per\n")
        f.write("column approximately match the real cohort where the real file\n")
        f.write("was available at generation time; joint dependence between\n")
        f.write("variables is intentionally not preserved.\n\n")
        f.write("Suitable for: pipeline smoke testing, schema validation,\n")
        f.write("adaptation to other cohorts.\n")
        f.write("NOT suitable for: drawing scientific conclusions.\n")
    print(f"  wrote {PROV}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
