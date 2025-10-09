#!/usr/bin/env python3
"""
Analyze SARS-CoV-2 aligned genomes vs Wuhan reference (NC_045512.2).

Features
- Per-sequence QC metrics: ambiguous bases, SNP/INDEL counts (+ Spike only)
- SNP list (long format)
- Sparse SNP presence/absence matrix (samples × sites)
- PCA + k-means clustering
  * k can be fixed (e.g., 3) or "auto" (silhouette-based selection over [kmin, kmax])
  * When k="auto", a k_selection.tsv is written with scores for each k
- Spike amino-acid mutation list

Inputs
- --aln  : aligned multi-FASTA (all sequences aligned to/ref length; e.g., Nextalign output)
- --ref  : reference FASTA (NC_045512.2) single-record; length must match alignment length
- --out  : output directory

Optional
- --k            : integer (e.g., 3) or "auto" (default: auto)
- --kmin/--kmax  : k search range when k="auto" (default: 2..10)
- --pca-max-dims : upper bound for PCA dimensions (default: 10)

Outputs (written to --out)
- stats_per_sequence.tsv
- mutations_long.tsv
- mut_matrix_sparse.npz               (samples × SNP sites; CSR)
- mut_matrix_sites.tsv                (SNP genomic sites)
- mut_matrix_samples.tsv              (sample IDs in matrix order)
- pca_kmeans.tsv                      (PCs + chosen k + cluster labels)
- spike_aa_mutations.tsv
- k_selection.tsv                     (only when k="auto")

Notes
- Insertions relative to reference are typically not represented in aligned FASTA columns
  from Nextalign/Nextclade; this script focuses on SNPs + deletions in the alignment.
"""

import sys
import os
import argparse
from collections import defaultdict
from typing import List, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.sparse import csr_matrix, save_npz
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from tqdm import tqdm

# Wuhan reference genomic coords for Spike (NC_045512.2)
SPIKE_START = 21563  # 1-based inclusive
SPIKE_END   = 25384  # 1-based inclusive


# ----------------------------- IO helpers -----------------------------

def read_first_record(fasta_path: str):
    for rec in SeqIO.parse(fasta_path, "fasta"):
        return rec
    raise RuntimeError(f"No sequences in {fasta_path}")

def read_all(fasta_path: str):
    return list(SeqIO.parse(fasta_path, "fasta"))


# ----------------------------- mutation logic -----------------------------

def count_ambiguous(seq_str: str) -> int:
    """Count IUPAC ambiguous bases (incl. N) and gaps as ambiguous (QC metric)."""
    amb = set("RYKMSWBDHVNrykmswbdhvn")
    return sum(1 for b in seq_str if (b in amb or b == '-'))

def within_spike(pos_1based: int) -> bool:
    return SPIKE_START <= pos_1based <= SPIKE_END

def find_mutations_vs_ref(ref_aln: str, query_aln: str) -> List[dict]:
    """
    Compare aligned reference string (same length) to aligned query.
    Reports SNPs and contiguous deletions relative to reference. Insertions are
    not in Nextalign/Nextclade aligned columns and are ignored here.
    """
    assert len(ref_aln) == len(query_aln)
    muts = []
    rpos = 0  # 1-based reference coordinate (increments when ref is not a gap)
    i = 0
    L = len(ref_aln)
    while i < L:
        r = ref_aln[i]
        q = query_aln[i]

        if r != '-':
            rpos += 1

        # identical / both gaps / ambiguous N in either -> skip
        if r == q or r.upper() == 'N' or q.upper() == 'N' or (r == '-' and q == '-'):
            i += 1
            continue

        # SNP
        if r != '-' and q != '-' and r != q:
            muts.append({"type": "SNP", "pos": rpos, "ref": r, "alt": q})
            i += 1
            continue

        # Deletion relative to reference (contiguous run of q == '-' while ref != '-')
        if r != '-' and q == '-':
            j = i
            delref = []
            start_pos = rpos
            while j < L and ref_aln[j] != '-' and query_aln[j] == '-':
                delref.append(ref_aln[j])
                j += 1
            muts.append({"type": "DEL", "pos": start_pos, "ref": "".join(delref), "alt": "-"})
            # rpos already incremented at loop head; add remaining consumed ref length - 1
            rpos += (len(delref) - 1)
            i = j
            continue

        # Insertions (q != '-' and r == '-') would require gapped ref; not typical here.
        i += 1
    return muts

def translate_spike(ref_aln: str, query_aln: str) -> Tuple[str, str]:
    """
    Extract and translate Spike coding region from aligned sequences.
    Keep only columns where BOTH ref and query are not gaps and within Spike coords.
    """
    ref_spike_nt = []
    qry_spike_nt = []
    ref_pos = 0  # 1-based genomic coordinate
    for r, q in zip(ref_aln, query_aln):
        if r != '-':
            ref_pos += 1
        if within_spike(ref_pos) and r != '-' and q != '-':
            ref_spike_nt.append(r)
            qry_spike_nt.append(q)
    ref_spike_seq = Seq("".join(ref_spike_nt))
    qry_spike_seq = Seq("".join(qry_spike_nt))
    return str(ref_spike_seq.translate(to_stop=False)), str(qry_spike_seq.translate(to_stop=False))

def aa_changes(ref_aa: str, qry_aa: str) -> List[str]:
    changes = []
    for i, (a, b) in enumerate(zip(ref_aa, qry_aa), start=1):
        if a != b:
            changes.append(f"{a}{i}{b}")
    return changes


# ----------------------------- PCA / k-means -----------------------------

def choose_k_auto(X: np.ndarray, kmin: int, kmax: int, random_state: int = 0):
    """
    Pick k in [kmin, kmax] by maximizing silhouette score on PCA-reduced data (up to 5 dims).
    Returns: chosen_k, results_df (k, silhouette, inertia).
    """
    results = []
    # Use min(5, n_features, n_samples-1) PCA dims to denoise
    n_samples, n_features = X.shape
    ncomp = max(1, min(5, n_features, n_samples - 1))

    # Robust PCA in case of tiny data
    try:
        Z = PCA(n_components=ncomp, random_state=random_state).fit_transform(X)
    except Exception:
        Z = X.copy()

    for k in range(kmin, kmax + 1):
        if k <= 1 or k >= n_samples:
            continue
        km = KMeans(n_clusters=k, n_init=10, random_state=random_state).fit(Z)
        labels = km.labels_
        # Silhouette requires >1 cluster and no singleton for stable score; catch failures
        try:
            sil = silhouette_score(Z, labels, metric="euclidean")
        except Exception:
            sil = np.nan
        inertia = float(km.inertia_)
        results.append({"k": k, "silhouette": sil, "inertia": inertia})

    if not results:
        # Fallback: no valid k tried
        return None, pd.DataFrame(columns=["k", "silhouette", "inertia"])

    df = pd.DataFrame(results)
    # Choose by highest silhouette; if ties/NaNs, use smallest inertia among the top silhouette rows
    # Drop NaN silhouettes for ranking, but keep them in the report
    valid = df.dropna(subset=["silhouette"])
    if len(valid) == 0:
        # If all NaN, choose k with minimal inertia
        best_row = df.loc[df["inertia"].idxmin()]
    else:
        max_sil = valid["silhouette"].max()
        candidates = valid.loc[valid["silhouette"] == max_sil]
        best_row = candidates.loc[candidates["inertia"].idxmin()]
    chosen_k = int(best_row["k"])
    return chosen_k, df


# ----------------------------- main -----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln", required=True, help="aligned multi-FASTA (Nextalign/Nextclade; same length as ref)")
    ap.add_argument("--ref", required=True, help="reference FASTA (Wuhan NC_045512.2), single record")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--k", default="auto", help='k for k-means (integer or "auto")')
    ap.add_argument("--kmin", type=int, default=2, help="min k when k='auto'")
    ap.add_argument("--kmax", type=int, default=10, help="max k when k='auto'")
    ap.add_argument("--pca-max-dims", type=int, default=10, help="upper bound for PCA components")
    ap.add_argument("--random-state", type=int, default=0, help="random seed for PCA/KMeans")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # Load reference and alignment
    ref_rec = read_first_record(args.ref)
    ref_nt = str(ref_rec.seq).upper().replace("\n", "")
    aln_records = read_all(args.aln)
    if not aln_records:
        sys.exit("No sequences in aligned FASTA.")

    aln_len = len(aln_records[0].seq)
    if any(len(rec.seq) != aln_len for rec in aln_records):
        sys.exit("Aligned sequences must all have the same length. Use Nextalign/Nextclade or MAFFT --keeplength.")

    if len(ref_nt) != aln_len:
        sys.exit(
            f"Reference length ({len(ref_nt)}) != aligned sequence length ({aln_len}). "
            "Ensure you aligned to/padded to reference coordinates."
        )

    ref_aln = ref_nt  # reference is ungapped but same length as alignment columns

    samples = [(rec.id, str(rec.seq).upper()) for rec in aln_records]

    # ---------- Per-sample stats and mutation listing ----------
    rows_stats = []
    mut_rows = []
    snp_sites = set()

    for sid, seq in tqdm(samples, desc="Scanning mutations"):
        n_amb = count_ambiguous(seq)
        muts = find_mutations_vs_ref(ref_aln, seq)

        n_snp = sum(1 for m in muts if m["type"] == "SNP")
        # INS are not produced by the alignment routine we assume; keep count for completeness
        n_ins = sum(1 for m in muts if m["type"] == "INS")
        n_del = sum(1 for m in muts if m["type"] == "DEL")
        n_snp_spike = sum(1 for m in muts if m["type"] == "SNP" and within_spike(m["pos"]))
        n_indel_spike = sum(1 for m in muts if m["type"] in ("INS", "DEL") and within_spike(m["pos"]))

        rows_stats.append({
            "sample": sid,
            "ambiguous_bases": n_amb,
            "snps": n_snp,
            "insertions": n_ins,
            "deletions": n_del,
            "snps_spike": n_snp_spike,
            "indels_spike": n_indel_spike
        })
        for m in muts:
            mut_rows.append({"sample": sid, **m})
            if m["type"] == "SNP":
                snp_sites.add(m["pos"])

    pd.DataFrame(rows_stats).to_csv(os.path.join(args.out, "stats_per_sequence.tsv"), sep="\t", index=False)
    pd.DataFrame(mut_rows).to_csv(os.path.join(args.out, "mutations_long.tsv"), sep="\t", index=False)

    # ---------- Sparse SNP matrix ----------
    snp_sites = sorted(snp_sites)
    site_index = {pos: i for i, pos in enumerate(snp_sites)}
    indptr = [0]; indices = []; data = []; sample_ids = []

    mut_df = pd.DataFrame(mut_rows)
    mut_df = mut_df[mut_df["type"] == "SNP"].copy()

    for sid, _ in samples:
        sample_ids.append(sid)
        s_sites = mut_df.loc[mut_df["sample"] == sid, "pos"].tolist()
        s_cols = sorted(site_index[p] for p in set(s_sites) if p in site_index)
        indices.extend(s_cols)
        data.extend([1] * len(s_cols))
        indptr.append(len(indices))

    X = csr_matrix((data, indices, indptr), shape=(len(sample_ids), len(snp_sites)), dtype=np.uint8)
    save_npz(os.path.join(args.out, "mut_matrix_sparse.npz"), X)
    pd.DataFrame({"site": snp_sites}).to_csv(os.path.join(args.out, "mut_matrix_sites.tsv"), sep="\t", index=False)
    pd.DataFrame({"sample": sample_ids}).to_csv(os.path.join(args.out, "mut_matrix_samples.tsv"), sep="\t", index=False)

    # ---------- PCA + k-means ----------
    # Guard: not enough features or samples
    if X.shape[1] == 0 or X.shape[0] < 2:
        # Write minimal pca_kmeans.tsv with no PCs and cluster=-1
        out_df = pd.DataFrame({"sample": sample_ids})
        out_df["kmeans_k"] = -1
        out_df["cluster"] = -1
        out_df.to_csv(os.path.join(args.out, "pca_kmeans.tsv"), sep="\t", index=False)
        # Nothing to choose for k
        # Still compute Spike AA mutations below
    else:
        X_dense = X.toarray().astype(float)
        # Determine PCA dimension conservatively
        ncomp = max(1, min(args.pca_max_dims, min(X_dense.shape) - 1))
        pca = PCA(n_components=ncomp, random_state=args.random_state)
        Z = pca.fit_transform(X_dense)

        # Choose k
        chosen_k = None
        k_selection_df = None

        if isinstance(args.k, str) and args.k.lower() == "auto":
            chosen_k, k_selection_df = choose_k_auto(Z, args.kmin, args.kmax, random_state=args.random_state)
            # If auto fails (e.g., too few samples), fall back to k=2 if possible
            if chosen_k is None:
                chosen_k = 2 if X.shape[0] > 2 else 1
        else:
            try:
                chosen_k = int(args.k)
            except ValueError:
                sys.exit("--k must be an integer or 'auto'.")

        # Train final k-means on top PCs (use first up to 5 dims or available)
        use_dims = min(5, Z.shape[1])
        if chosen_k <= 1 or chosen_k > max(1, Z.shape[0] - 1):
            # Degenerate case: cannot cluster
            labels = np.full(Z.shape[0], -1, dtype=int)
            chosen_k_out = -1
        else:
            km = KMeans(n_clusters=chosen_k, n_init=10, random_state=args.random_state)
            labels = km.fit_predict(Z[:, :use_dims])
            chosen_k_out = chosen_k

        out_df = pd.DataFrame({"sample": sample_ids})
        for i in range(Z.shape[1]):
            out_df[f"PC{i+1}"] = Z[:, i]
        out_df["kmeans_k"] = chosen_k_out
        out_df["cluster"] = labels
        out_df.to_csv(os.path.join(args.out, "pca_kmeans.tsv"), sep="\t", index=False)

        # Write k_selection.tsv only when k="auto" and we computed scores
        if (isinstance(args.k, str) and args.k.lower() == "auto") and (k_selection_df is not None):
            k_selection_df.to_csv(os.path.join(args.out, "k_selection.tsv"), sep="\t", index=False)

    # ---------- Spike amino-acid mutations per sample ----------
    aa_rows = []
    for sid, seq in tqdm(samples, desc="Translating Spike"):
        ref_aa, qry_aa = translate_spike(ref_aln, seq)
        for aamut in aa_changes(ref_aa, qry_aa):
            aa_rows.append({"sample": sid, "Spike_AA_mut": aamut})
    pd.DataFrame(aa_rows).to_csv(os.path.join(args.out, "spike_aa_mutations.tsv"), sep="\t", index=False)

    print("\nDone. Outputs written to:", args.out)


if __name__ == "__main__":
    main()


