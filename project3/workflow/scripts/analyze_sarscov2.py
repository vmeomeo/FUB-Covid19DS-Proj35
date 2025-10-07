#!/usr/bin/env python3
import sys, os, math
from collections import defaultdict, Counter
from typing import List, Tuple, Dict

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.sparse import csr_matrix, save_npz
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from tqdm import tqdm
import argparse

"""
USAGE:
  python analyze_sarscov2.py \
    --aln work/aligned.fasta \
    --ref data/NC_045512.2.fasta \
    --out outdir

Inputs:
- aligned.fasta : Nextclade/Nextalign/MAFFT output. All sequences must have SAME length as reference.
- ref fasta     : Wuhan reference (NC_045512.2) single-record FASTA

Outputs written to --out:
- stats_per_sequence.tsv
- mutations_long.tsv
- mut_matrix_sparse.npz  (samples x variant-sites, 0/1) + helper TSVs
- pca_kmeans.tsv
- spike_aa_mutations.tsv
"""

SPIKE_START = 21563  # 1-based genomic coords on NC_045512.2
SPIKE_END   = 25384  # inclusive

def read_first_record(fasta):
    for rec in SeqIO.parse(fasta, "fasta"):
        return rec
    raise RuntimeError(f"No sequences in {fasta}")

def read_all(fasta):
    return list(SeqIO.parse(fasta, "fasta"))

def count_ambiguous(seq_str: str) -> int:
    # Count IUPAC ambiguity (including N) and treat gaps as ambiguous for a simple QC metric
    amb = set("RYKMSWBDHVNn")
    return sum(1 for b in seq_str if (b in amb or b == '-'))

def within_spike(pos_1based:int) -> bool:
    return SPIKE_START <= pos_1based <= SPIKE_END

def find_mutations_vs_ref(ref_aln: str, query_aln: str) -> List[dict]:
    """
    Compare aligned reference string (same length as query_aln) to aligned query.
    Reports SNPs and small deletions relative to ref; insertions are not present
    in Nextclade/Nextalign aligned FASTA and are typically provided separately.
    """
    assert len(ref_aln) == len(query_aln)
    muts = []
    rpos = 0  # 1-based reference coordinate (increments when ref is not a gap)
    i = 0
    while i < len(ref_aln):
        r = ref_aln[i]
        q = query_aln[i]

        if r != '-':
            rpos += 1

        # skip perfect match or N/ambiguous comparisons
        if r == q or r.upper() == 'N' or q.upper() == 'N' or (r == '-' and q == '-'):
            i += 1
            continue

        # SNP
        if r != '-' and q != '-' and r != q:
            muts.append({"type":"SNP","pos":rpos,"ref":r,"alt":q})
            i += 1
            continue

        # Deletion relative to reference (contiguous run)
        if r != '-' and q == '-':
            j = i
            delref = []
            start_pos = rpos
            # consume a contiguous deletion run
            while j < len(ref_aln) and ref_aln[j] != '-' and query_aln[j] == '-':
                delref.append(ref_aln[j])
                # rpos will be incremented at loop head; account after loop
                j += 1
            muts.append({"type":"DEL","pos":start_pos, "ref":"".join(delref), "alt":"-"})
            # advance rpos by remaining consumed ref letters minus the auto-increment already done
            rpos += (len(delref) - 1)
            i = j
            continue

        # Insertions relative to reference would need a gapped ref ('-') column.
        # Nextclade/Nextalign do not put insertions into the aligned FASTA columns,
        # so we ignore insertions here. They are in nextclade's insertions.tsv if needed.

        i += 1
    return muts

def translate_spike(ref_aln: str, query_aln: str) -> Tuple[str, str]:
    """
    Translate Spike coding region from aligned sequences.
    Keep only columns where BOTH ref and query are not gaps and within Spike.
    """
    ref_spike_nt = []
    qry_spike_nt = []
    ref_pos = 0  # 1-based genomic position
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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln", required=True, help="aligned multi-FASTA (all seqs same length as reference)")
    ap.add_argument("--ref", required=True, help="reference FASTA (NC_045512.2) single-record")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--k", type=int, default=6, help="k-means clusters")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # Read reference (ungapped) and alignment (gapped to ref length)
    ref_nt = str(read_first_record(args.ref).seq).upper().replace("\n","")
    aln_records = read_all(args.aln)
    if not aln_records:
        sys.exit("No sequences in aligned FASTA.")
    aln_len = len(aln_records[0].seq)
    if any(len(rec.seq) != aln_len for rec in aln_records):
        sys.exit("Aligned sequences have different lengths. Ensure you used Nextclade/Nextalign or MAFFT --keeplength.")
    if len(ref_nt) != aln_len:
        # With Nextclade/Nextalign, aligned length should equal reference length
        sys.exit(f"Reference length ({len(ref_nt)}) != aligned sequence length ({aln_len}). "
                 "Use Nextclade/Nextalign or pad to reference coordinates.")

    ref_aln = ref_nt  # same length as alignment columns (no gaps in ref)

    # Build list of samples from alignment (use all records; ref is NOT present here)
    samples = [(rec.id, str(rec.seq).upper()) for rec in aln_records]

    # Per-sample stats and mutation listing
    rows_stats = []
    mut_rows = []
    snp_sites = set()

    for sid, seq in tqdm(samples, desc="Scanning mutations"):
        n_amb = count_ambiguous(seq)
        muts = find_mutations_vs_ref(ref_aln, seq)
        n_snp = sum(1 for m in muts if m["type"] == "SNP")
        n_ins = sum(1 for m in muts if m["type"] == "INS")
        n_del = sum(1 for m in muts if m["type"] == "DEL")
        n_snp_spike = sum(1 for m in muts if m["type"] == "SNP" and within_spike(m["pos"]))
        n_indel_spike = sum(1 for m in muts if m["type"] in ("INS","DEL") and within_spike(m["pos"]))

        rows_stats.append({
            "sample": sid,
            "ambiguous_bases": n_amb,
            "snps": n_snp, "insertions": n_ins, "deletions": n_del,
            "snps_spike": n_snp_spike, "indels_spike": n_indel_spike
        })
        for m in muts:
            mut_rows.append({"sample": sid, **m})
            if m["type"] == "SNP":
                snp_sites.add(m["pos"])

    pd.DataFrame(rows_stats).to_csv(os.path.join(args.out, "stats_per_sequence.tsv"), sep="\t", index=False)
    pd.DataFrame(mut_rows).to_csv(os.path.join(args.out, "mutations_long.tsv"), sep="\t", index=False)

    # Sparse SNP matrix (samples x SNP sites)
    snp_sites = sorted(snp_sites)
    site_index = {pos:i for i, pos in enumerate(snp_sites)}
    indptr = [0]; indices = []; data = []; sample_ids = []

    mut_df = pd.DataFrame(mut_rows)
    mut_df = mut_df[mut_df["type"] == "SNP"].copy()

    for sid, _ in samples:
        sample_ids.append(sid)
        s_sites = mut_df.loc[mut_df["sample"] == sid, "pos"].tolist()  # <-- fixed column access
        s_cols = sorted(site_index[p] for p in set(s_sites) if p in site_index)
        indices.extend(s_cols)
        data.extend([1] * len(s_cols))
        indptr.append(len(indices))

    X = csr_matrix((data, indices, indptr), shape=(len(sample_ids), len(snp_sites)), dtype=np.uint8)
    save_npz(os.path.join(args.out, "mut_matrix_sparse.npz"), X)
    pd.DataFrame({"site": snp_sites}).to_csv(os.path.join(args.out, "mut_matrix_sites.tsv"), sep="\t", index=False)
    pd.DataFrame({"sample": sample_ids}).to_csv(os.path.join(args.out, "mut_matrix_samples.tsv"), sep="\t", index=False)

    # PCA + k-means (guard for small problems)
    # If there are no SNP sites or fewer than 2 samples, skip PCA/kmeans gracefully.
    if X.shape[1] == 0 or X.shape[0] < 2:
        out_df = pd.DataFrame({"sample": sample_ids})
        out_df["kmeans_k"] = args.k
        out_df["cluster"] = -1
        out_df.to_csv(os.path.join(args.out, "pca_kmeans.tsv"), sep="\t", index=False)
    else:
        X_dense = X.toarray().astype(float)
        ncomp = max(1, min(10, min(X_dense.shape) - 1))
        pca = PCA(n_components=ncomp, random_state=0)
        Z = pca.fit_transform(X_dense)
        use_dims = min(5, ncomp)
        # k must be <= number of samples; cap it
        k = min(args.k, max(1, X.shape[0] - 1))
        km = KMeans(n_clusters=k, n_init="auto", random_state=0).fit(Z[:, :use_dims])

        out_df = pd.DataFrame({"sample": sample_ids})
        for i in range(ncomp):
            out_df[f"PC{i+1}"] = Z[:, i]
        out_df["kmeans_k"] = k
        out_df["cluster"] = km.labels_
        out_df.to_csv(os.path.join(args.out, "pca_kmeans.tsv"), sep="\t", index=False)

    # Spike amino-acid mutations per sample
    aa_rows = []
    for sid, seq in tqdm(samples, desc="Translating Spike"):
        ref_aa, qry_aa = translate_spike(ref_aln, seq)
        for aamut in aa_changes(ref_aa, qry_aa):
            aa_rows.append({"sample": sid, "Spike_AA_mut": aamut})
    pd.DataFrame(aa_rows).to_csv(os.path.join(args.out, "spike_aa_mutations.tsv"), sep="\t", index=False)

    print("\nDone. Outputs in", args.out)

if __name__ == "__main__":
    main()

