
# Project3                            

This repository contains a Snakemake pipeline created by Kompal Fayyaz and Anh Vo.

-fetch genomes,
-align them to Wuhan-Hu-1 (NC_045512.2),
-call differences (SNPs/INDELs) vs reference,
-extract Spike mutations & simple per-sample stats,
-(optional) merge Bloom-lab escape scores,
-build a mutation matrix, PCA + k-means clusters, and
-generate quick figures.

It is designed to run on a laptop or HPC with conda/mamba.

# Features

-Robust input handling for FASTA consensus genomes (not FASTQ reads).
-Choice of aligner: Nextclade/Nextalign (default) or MAFFT (chunked, low-RAM).
-Clean, interpretable outputs (TSV + PNG).
-Reproducible via conda envs pinned in envs/.

# Repository layout
project3/
├── data/
│   ├── mut_profiles                 # GISAID mutation profiles
│   └── ref/NC_045512.2.fasta        # auto-downloaded by rule
│
├── workflow/
│   ├── envs/
│   │   ├── base.yaml
│   │   └── minimap2.yaml
│   ├── rules/
│   │   ├── sth1.smk
│   │   └── sth2.smk
│   ├── scripts/
│   │   ├── gisaid_download.R
│   │   ├── parse_minimap2_cs.py
│   │   └── pca_kmeans_and_risk.py
│   └── Snakefile
│
└── out/                             # results land here

project3/
├─ Snakefile

├─ README.md

├─ guide_proj3.MD

├─ config/

│  └─ config.yaml

├─ data/  # place your input FASTA(s) here

├─ out/   # results go here (can change in config)

├─ workflow/

│  ├─ rules/

│  │  ├─ get_reference.smk   # fetch NC_045512.2 if missing

│  │  ├─ io.smk              # sanity checks, paths

│  │  ├─ align_nextalign.smk # nextclade/nextalign alignment (default)

│  │  ├─ align.smk           # MAFFT (chunked) alternative

│  │  ├─ analyze.smk         # run analyzers (TSVs, PCA/kmeans)

│  │  └─ bloom.smk           # optional: Bloom escape integration

│  └─ envs/
│     ├─ base.yaml           # python: pandas, biopython, sklearn, matplotlib, tqdm

│     ├─ nextclade.yaml      # nextclade / nextalign CLI

│     └─ mafft.yaml          # mafft, seqkit
└─ scripts/
   ├─ analyze_sarscov2.py    # differences, stats, mutation matrix, PCA/kmeans, Spike AAs
   
   ├─ visualize_sarscov2.py  # quick PNG figures
   
   └─ compute_bloom_risk.py  # (optional) join Bloom escape map


# Inputs

Consensus genomes (FASTA): path defined in config.yaml (input_fasta), e.g. data/risk-assessment-sequences.fasta

Reference (Wuhan-Hu-1): fetched automatically (NC_045512.2) or provide local path via reference_fasta


# Outputs

All results land under results/:

results/aligned.fasta — aligned genomes to reference length

results/insertions.tsv — (nextclade/nextalign) insertions summary (if available)

results/stats_per_sequence.tsv — Ns, SNPs, INDELs (+ Spike-only) per genome

results/mutations_long.tsv — long table of all diffs vs reference

results/mut_matrix_sparse.npz + mut_matrix_sites.tsv + mut_matrix_samples.tsv

results/pca_kmeans.tsv — PCA components + cluster labels

results/spike_aa_mutations.tsv — AA changes in Spike (e.g., E484K)

results/stats_with_bloom.tsv — (optional) adds simple Bloom escape burden

figures/*.png — histograms, PCA scatter, per-site SNP frequency plots
#


# Reproducibility

All software is declared in envs/*.yaml.
All parameters are in config.yaml.


The DAG (snakemake --dag) documents the exact dependency graph.
Set a fixed kmeans_k and seeds (already fixed in the Python script) for stable clusters.

# Citations         

Nextclade / Nextalign (Nextstrain)
MAFFT alignment
Pandas / scikit-learn / Biopython / matplotlib
Bloom Lab (optional escape maps)
NCBI (Wuhan reference NC_045512.2)
Please cite these tools if you publish results.



License
Add your preferred license (e.g., MIT) in LICENSE.

