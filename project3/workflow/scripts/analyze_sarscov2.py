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
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from tqdm import tqdm
import argparse
import matplotlib.pyplot as plt

SPIKE_START = 21563
SPIKE_END   = 25384

def read_first_record(fasta):
    for rec in SeqIO.parse(fasta, "fasta"):
        return rec
    raise RuntimeError(f"No sequences in {fasta}")

def read_all(fasta):
    return list(SeqIO.parse(fasta, "fasta"))

def count_ambiguous(seq_str: str) -> int:
    amb = set("RYKMSWBDHVNn")
    return sum(1 for b in seq_str if (b in amb or b == '-'))

def within_spike(pos_1based:int) -> bool:
    return SPIKE_START <= pos_1based <= SPIKE_END

def find_mutations_vs_ref(ref_aln: str, query_aln: str) -> List[dict]:
    assert len(ref_aln) == len(query_aln)
    muts = []
    rpos = 0
    i = 0
    while i < len(ref_aln):
        r = ref_aln[i]; q = query_aln[i]
        if r != '-':
            rpos += 1
        if r == q or r.upper() == 'N' or q.upper() == 'N' or (r == '-' and q == '-'):
            i += 1; continue
        if r != '-' and q != '-' and r != q:
            muts.append({"type":"SNP","pos":rpos,"ref":r,"alt":q})
            i += 1; continue
        if r != '-' and q == '-':
            j = i; delref = []; start_pos = rpos
            while j < len(ref_aln) and ref_aln[j] != '-' and query_aln[j] == '-':
                delref.append(ref_aln[j]); j += 1
            muts.append({"type":"DEL","pos":start_pos, "ref":"".join(delref), "alt":"-"})
            rpos += (len(delref) - 1); i = j; continue
        i += 1
    return muts

def translate_spike(ref_aln: str, query_aln: str) -> Tuple[str, str]:
    ref_spike_nt = []; qry_spike_nt = []; ref_pos = 0
    for r, q in zip(ref_aln, query_aln):
        if r != '-': ref_pos += 1
        if within_spike(ref_pos) and r != '-' and q != '-':
            ref_spike_nt.append(r); qry_spike_nt.append(q)
    ref_spike_seq = Seq("".join(ref_spike_nt))
    qry_spike_seq = Seq("".join(qry_spike_nt))
    return str(ref_spike_seq.translate(to_stop=False)), str(qry_spike_seq.translate(to_stop=False))

def aa_changes(ref_aa: str, qry_aa: str) -> List[str]:
    return [f"{a}{i}{b}" for i,(a,b) in enumerate(zip(ref_aa, qry_aa), start=1) if a != b]

def choose_k_and_fit(Z, k_arg, kmin=2, kmax=10, dims=5, select="silhouette", random_state=0):
    """
    If k_arg is an int >=2: use that fixed k.
    If k_arg is None or 'auto' or <=1: auto-select k in [kmin, kmax] by selected metric.
    Returns (labels, chosen_k, metrics_df, km_model)
    """
    X = Z[:, :dims]
    n = X.shape[0]
    kmax = min(kmax, max(2, n-1))
    kmin = min(max(2, kmin), kmax)

    def fit_k(k):
        km = KMeans(n_clusters=k, n_init="auto", random_state=random_state).fit(X)
        lab = km.labels_
        # metrics
        try: sil = silhouette_score(X, lab)
        except Exception: sil = np.nan
        ch = calinski_harabasz_score(X, lab)
        db = davies_bouldin_score(X, lab)
        return km, lab, sil, ch, db

    # Fixed k?
    if isinstance(k_arg, int) and k_arg >= 2:
        k = min(k_arg, kmax)
        km, lab, sil, ch, db = fit_k(k)
        m = pd.DataFrame([{"k":k, "silhouette":sil, "calinski_harabasz":ch, "davies_bouldin":db, "inertia":km.inertia_}])
        return lab, k, m, km

    # Auto mode
    rows = []; best = None
    for k in range(kmin, kmax+1):
        km, lab, sil, ch, db = fit_k(k)
        rows.append({"k":k, "silhouette":sil, "calinski_harabasz":ch, "davies_bouldin":db, "inertia":km.inertia_})
        # choose by selected metric, with sensible tiebreaks
        if select == "ch":
            key = (ch, np.nan_to_num(sil, nan=-1), -db)
        elif select == "db":
            key = (-db, np.nan_to_num(sil, nan=-1), ch)
        else:  # silhouette default
            key = (np.nan_to_num(sil, nan=-1), ch, -db)
        if (best is None) or (key > best[0]):
            best = (key, k, km, lab)
    metrics = pd.DataFrame(rows)
    _, k_star, km_star, lab_star = best
    return lab_star, k_star, metrics, km_star

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", default="auto", help="int (fixed k) or 'auto' (choose k)")
    ap.add_argument("--kmin", type=int, default=2)
    ap.add_argument("--kmax", type=int, default=10)
    ap.add_argument("--k_select", choices=["silhouette","ch","db"], default="silhouette",
                    help="metric to pick k when --k auto")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # Read and sanity checks
    ref_nt = str(read_first_record(args.ref).seq).upper().replace("\n","")
    aln_records = read_all(args.aln)
    if not aln_records: sys.exit("No sequences in aligned FASTA.")
    aln_len = len(aln_records[0].seq)
    if any(len(rec.seq) != aln_len for rec in aln_records):
        sys.exit("Aligned sequences have different lengths; use Nextalign/Nextclade or MAFFT --keeplength.")
    if len(ref_nt) != aln_len:
        sys.exit(f"Reference length ({len(ref_nt)}) != aligned length ({aln_len}).")

    ref_aln = ref_nt
    samples = [(rec.id, str(rec.seq).upper()) for rec in aln_records]

    # Per-sample stats + long mutations + SNP site set
    rows_stats=[]; mut_rows=[]; snp_sites=set()
    for sid, seq in tqdm(samples, desc="Scanning mutations"):
        n_amb = count_ambiguous(seq)
        muts = find_mutations_vs_ref(ref_aln, seq)
        n_snp = sum(1 for m in muts if m["type"]=="SNP")
        n_ins = sum(1 for m in muts if m["type"]=="INS")
        n_del = sum(1 for m in muts if m["type"]=="DEL")
        n_snp_spike = sum(1 for m in muts if m["type"]=="SNP" and within_spike(m["pos"]))
        n_indel_spike = sum(1 for m in muts if m["type"] in ("INS","DEL") and within_spike(m["pos"]))
        rows_stats.append({"sample":sid,"ambiguous_bases":n_amb,"snps":n_snp,"insertions":n_ins,
                           "deletions":n_del,"snps_spike":n_snp_spike,"indels_spike":n_indel_spike})
        for m in muts:
            mut_rows.append({"sample":sid, **m})
            if m["type"]=="SNP": snp_sites.add(m["pos"])

    pd.DataFrame(rows_stats).to_csv(os.path.join(args.out,"stats_per_sequence.tsv"), sep="\t", index=False)
    pd.DataFrame(mut_rows).to_csv(os.path.join(args.out,"mutations_long.tsv"), sep="\t", index=False)

    # Sparse SNP matrix
    snp_sites = sorted(snp_sites)
    site_index = {p:i for i,p in enumerate(snp_sites)}
    indptr=[0]; indices=[]; data=[]; sample_ids=[]

    mut_df = pd.DataFrame(mut_rows)
    mut_df = mut_df[mut_df["type"]=="SNP"].copy()

    for sid,_ in samples:
        sample_ids.append(sid)
        s_sites = mut_df.loc[mut_df["sample"]==sid, "pos"].tolist()
        s_cols = sorted(site_index[p] for p in set(s_sites) if p in site_index)
        indices.extend(s_cols); data.extend([1]*len(s_cols)); indptr.append(len(indices))

    X = csr_matrix((data, indices, indptr), shape=(len(sample_ids), len(snp_sites)), dtype=np.uint8)
    save_npz(os.path.join(args.out,"mut_matrix_sparse.npz"), X)
    pd.DataFrame({"site":snp_sites}).to_csv(os.path.join(args.out,"mut_matrix_sites.tsv"), sep="\t", index=False)
    pd.DataFrame({"sample":sample_ids}).to_csv(os.path.join(args.out,"mut_matrix_samples.tsv"), sep="\t", index=False)

    # PCA + KMeans (with auto-k option)
    if X.shape[1]==0 or X.shape[0]<2:
        out_df = pd.DataFrame({"sample":sample_ids})
        out_df["kmeans_k"] = -1
        out_df["cluster"] = -1
        out_df.to_csv(os.path.join(args.out,"pca_kmeans.tsv"), sep="\t", index=False)
    else:
        Xd = X.toarray().astype(float)
        ncomp = max(1, min(10, min(Xd.shape)-1))
        pca = PCA(n_components=ncomp, random_state=0)
        Z = pca.fit_transform(Xd)
        use_dims = min(5, ncomp)

        # parse k argument
        try:
            k_arg = int(args.k)
        except ValueError:
            k_arg = -1  # auto

        labels, k_star, metrics_df, km = choose_k_and_fit(
            Z, k_arg=k_arg, kmin=args.kmin, kmax=args.kmax,
            dims=use_dims, select=args.k_select, random_state=0
        )
        metrics_df.to_csv(os.path.join(args.out, "kmeans_metrics.tsv"), sep="\t", index=False)

        out_df = pd.DataFrame({"sample":sample_ids})
        for i in range(ncomp):
            out_df[f"PC{i+1}"] = Z[:, i]
        out_df["kmeans_k"] = k_star
        out_df["cluster"] = labels
        out_df.to_csv(os.path.join(args.out,"pca_kmeans.tsv"), sep="\t", index=False)

    # Spike AA mutations per sample
    aa_rows=[]
    for sid, seq in tqdm(samples, desc="Translating Spike"):
        ref_aa, qry_aa = translate_spike(ref_aln, seq)
        for aamut in aa_changes(ref_aa, qry_aa):
            aa_rows.append({"sample":sid, "Spike_AA_mut":aamut})
    pd.DataFrame(aa_rows).to_csv(os.path.join(args.out,"spike_aa_mutations.tsv"), sep="\t", index=False)

    print("\nDone. Outputs in", args.out)

if __name__ == "__main__":
    main()



