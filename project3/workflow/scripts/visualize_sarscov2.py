#!/usr/bin/env python3
import os, argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def load_tsv(path):
    return pd.read_csv(path, sep="\t")

def plot_mutation_load(stats, ofn):
    plt.figure(figsize=(7,4))
    vals = stats["snps"] + stats["insertions"] + stats["deletions"]
    plt.hist(vals, bins=20)
    plt.xlabel("Total mutations per genome (SNPs + INDELs)")
    plt.ylabel("Number of genomes")
    plt.title("Mutation load distribution")
    plt.tight_layout(); plt.savefig(ofn); plt.close()

def plot_ambiguous(stats, ofn):
    plt.figure(figsize=(7,4))
    plt.hist(stats["ambiguous_bases"], bins=20)
    plt.xlabel("Ambiguous bases (N/IUPAC + gaps in MSA)")
    plt.ylabel("Number of genomes")
    plt.title("Ambiguity distribution")
    plt.tight_layout(); plt.savefig(ofn); plt.close()

def plot_pca(pca_km, ofn):
    # handle the case where PCA was skipped
    if "PC1" not in pca_km.columns or "PC2" not in pca_km.columns:
        # make a tiny placeholder figure
        plt.figure(figsize=(6,4))
        plt.text(0.5,0.5,"PCA not available (no SNP sites or too few samples)", ha="center")
        plt.axis("off"); plt.tight_layout(); plt.savefig(ofn); plt.close(); return
    plt.figure(figsize=(6,5))
    clusters = pca_km["cluster"].astype(int)
    # a simple color per cluster:
    for c in sorted(clusters.unique()):
        mask = clusters==c
        plt.scatter(pca_km.loc[mask,"PC1"], pca_km.loc[mask,"PC2"], s=24, label=f"cluster {c}")
    plt.xlabel("PC1"); plt.ylabel("PC2"); plt.title("PCA of SNP matrix")
    plt.legend(loc="best", frameon=False)
    plt.tight_layout(); plt.savefig(ofn); plt.close()

def plot_top_spike_aas(spike_aa, ofn, topn=20):
    if spike_aa.empty:
        plt.figure(figsize=(6,4))
        plt.text(0.5,0.5,"No Spike AA mutations found", ha="center")
        plt.axis("off"); plt.tight_layout(); plt.savefig(ofn); plt.close(); return
    counts = spike_aa["Spike_AA_mut"].value_counts().head(topn)[::-1]
    plt.figure(figsize=(8, max(4, 0.3*len(counts))))
    plt.barh(counts.index, counts.values)
    plt.xlabel("Genomes with mutation"); plt.ylabel("Spike AA change")
    plt.title(f"Top {len(counts)} Spike amino-acid mutations")
    plt.tight_layout(); plt.savefig(ofn); plt.close()

def plot_snp_site_frequency(mut_long, ofn):
    # frequency across all genomes (SNPs only), by genomic position
    snps = mut_long[mut_long["type"]=="SNP"].copy()
    if snps.empty:
        plt.figure(figsize=(6,4))
        plt.text(0.5,0.5,"No SNPs found", ha="center"); plt.axis("off")
        plt.tight_layout(); plt.savefig(ofn); plt.close(); return
    freq = snps.groupby("pos").size().reset_index(name="count")
    plt.figure(figsize=(10,4))
    plt.scatter(freq["pos"], freq["count"], s=8)
    plt.xlabel("Genomic position (NC_045512.2)"); plt.ylabel("SNP count across genomes")
    plt.title("Per-site SNP frequency")
    plt.tight_layout(); plt.savefig(ofn); plt.close()

def plot_spike_only_snp_freq(mut_long, ofn, spike_start=21563, spike_end=25384):
    snps = mut_long[(mut_long["type"]=="SNP") &
                    (mut_long["pos"].between(spike_start, spike_end))].copy()
    if snps.empty:
        plt.figure(figsize=(6,4))
        plt.text(0.5,0.5,"No Spike SNPs found", ha="center"); plt.axis("off")
        plt.tight_layout(); plt.savefig(ofn); plt.close(); return
    freq = snps.groupby("pos").size().reset_index(name="count")
    plt.figure(figsize=(10,4))
    plt.scatter(freq["pos"], freq["count"], s=10)
    plt.xlabel("Spike genomic position"); plt.ylabel("SNP count across genomes")
    plt.title("Spike-region SNP frequency")
    plt.tight_layout(); plt.savefig(ofn); plt.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="indir", default="out", help="folder with TSV outputs")
    ap.add_argument("--figdir", default="figures", help="where to save PNGs")
    args = ap.parse_args()

    ensure_dir(args.figdir)

    stats = load_tsv(os.path.join(args.indir, "stats_per_sequence.tsv"))
    mut_long = load_tsv(os.path.join(args.indir, "mutations_long.tsv"))
    pca_km = load_tsv(os.path.join(args.indir, "pca_kmeans.tsv"))
    # spike AA file may be missing if something went wrong; handle gracefully
    spike_path = os.path.join(args.indir, "spike_aa_mutations.tsv")
    spike_aa = load_tsv(spike_path) if os.path.exists(spike_path) else pd.DataFrame(columns=["sample","Spike_AA_mut"])

    plot_mutation_load(stats, os.path.join(args.figdir, "hist_mutation_load.png"))
    plot_ambiguous(stats, os.path.join(args.figdir, "hist_ambiguous_bases.png"))
    plot_pca(pca_km, os.path.join(args.figdir, "pca_clusters.png"))
    plot_top_spike_aas(spike_aa, os.path.join(args.figdir, "top_spike_aas.png"))
    plot_snp_site_frequency(mut_long, os.path.join(args.figdir, "snp_site_frequency.png"))
    plot_spike_only_snp_freq(mut_long, os.path.join(args.figdir, "spike_snp_site_frequency.png"))

    print("Saved figures to:", args.figdir)

if __name__ == "__main__":
    main()

