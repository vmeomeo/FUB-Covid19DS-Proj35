#!/usr/bin/env python3
import os, io, math, sys
import pandas as pd
import numpy as np
import requests

# --- CONFIG (matches your earlier pipeline) ---
SPIKE_START = 21563          # NC_045512.2 coords (1-based)
SPIKE_END   = 25384
RBD_AA_START, RBD_AA_END = 319, 541
USE_BLOOM = True
BLOOM_RAW = "https://raw.githubusercontent.com/jbloomlab/SARS2-RBD-escape-calc/main/results/escape_chart_data.csv"

OUTDIR = "out"               # where stats/mutations live
STATS_TSV = os.path.join(OUTDIR, "stats_per_sequence.tsv")
MUTS_TSV  = os.path.join(OUTDIR, "mutations_long.tsv")
OUT_TSV   = os.path.join(OUTDIR, "stats_with_bloom.tsv")

def nt_to_spike_aa(nt_pos_1b: int) -> int:
    """Map genomic nt position -> 1-based amino-acid index in Spike."""
    return 1 + (nt_pos_1b - SPIKE_START) // 3

def load_bloom():
    if not USE_BLOOM:
        return None, None
    try:
        r = requests.get(BLOOM_RAW, timeout=60)
        r.raise_for_status()
        esc = pd.read_csv(io.BytesIO(r.content))
        # find a site and value column flexibly
        site_col = next((c for c in ["site","aa_site","position","aa_position"] if c in esc.columns), None)
        val_col  = next((c for c in ["escape","escape_score","median_escape"] if c in esc.columns), None)
        if site_col is None or val_col is None:
            print("[WARN] Bloom table loaded but site/value columns not recognized; ignoring Bloom weights.")
            return None, None
        print(f"[OK] Bloom escape table loaded. Using columns: site={site_col}, value={val_col}")
        # Pre-aggregate per AA site (mean)
        site_escape = esc.groupby(site_col)[val_col].mean().rename("bloom_escape").reset_index()
        site_escape = site_escape.rename(columns={site_col: "aa"})
        site_escape["aa"] = site_escape["aa"].astype(int)
        return site_escape, "bloom_escape"
    except Exception as e:
        print("[WARN] Bloom load failed:", e)
        return None, None

def main():
    # Load inputs
    if not (os.path.exists(STATS_TSV) and os.path.exists(MUTS_TSV)):
        sys.exit(f"Missing inputs. Expect {STATS_TSV} and {MUTS_TSV}")

    stats = pd.read_csv(STATS_TSV, sep="\t")
    muts  = pd.read_csv(MUTS_TSV, sep="\t")

    # keep SNPs in Spike coordinates
    snps = muts[(muts["type"] == "SNP") &
                (muts["pos"].between(SPIKE_START, SPIKE_END))].copy()

    if snps.empty:
        print("[INFO] No Spike SNPs found — writing passthrough with toy_risk= -0.1*ambiguous.")
        out = stats.copy()
        out["toy_risk"] = -0.1 * out["ambiguous_bases"].fillna(0)
        out.to_csv(OUT_TSV, sep="\t", index=False)
        print(f"[OK] Wrote {OUT_TSV}")
        print(out.sort_values("toy_risk", ascending=False).head(10))
        return

    # Map nt pos -> AA site; collapse to unique AA sites per sample (avoid double counting multiple SNPs in same codon)
    snps["aa"] = ((snps["pos"] - SPIKE_START) // 3 + 1).astype(int)
    per_sample_aa = (snps.groupby("sample")["aa"]
                          .apply(lambda s: sorted(set(int(x) for x in s)))
                          .reset_index(name="aa_sites"))

    # Load Bloom escape (optional)
    bloom_df, bloom_col = load_bloom()
    bloom_map = {}
    if bloom_df is not None:
        bloom_map = {int(row.aa): float(row[bloom_col]) for _, row in bloom_df.iterrows()}

    # Compute toy risk
    # base weight = 1 per AA site; + (Bloom escape at that site if available)
    # RBD AA sites get an extra *2 multiplier (or simple +1)? We'll follow your sketch: multiply by 2.0
    def aa_weight(aa_site: int) -> float:
        w = 1.0
        if RBD_AA_START <= aa_site <= RBD_AA_END:
            w *= 2.0
        if bloom_map:
            w += bloom_map.get(aa_site, 0.0)
        return w

    risk_rows = []
    ambig = stats.set_index("sample")["ambiguous_bases"].to_dict()

    for _, row in per_sample_aa.itertuples(index=False):
        pass  # placeholder to keep tuple unpacking doc alive

    for rec in per_sample_aa.itertuples(index=False):
        sid = rec.sample
        aa_sites = rec.aa_sites
        score = sum(aa_weight(a) for a in aa_sites)
        score -= 0.1 * float(ambig.get(sid, 0))  # tiny penalty for ambiguity
        risk_rows.append((sid, score))

    risk_df = pd.DataFrame(risk_rows, columns=["sample","toy_risk"])
    out = stats.merge(risk_df, on="sample", how="left")
    out["toy_risk"] = out["toy_risk"].fillna(-0.1 * out["ambiguous_bases"].fillna(0))

    out.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"[OK] Wrote {OUT_TSV}")
    print(out.sort_values("toy_risk", ascending=False).head(10))

if __name__ == "__main__":
    main()

