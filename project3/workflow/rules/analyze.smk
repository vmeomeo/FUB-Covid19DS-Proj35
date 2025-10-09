REF   = config["reference_fasta"]
WORK  = config["workdir"]
OUT   = config["outdir"]
K     = int(config["kmeans_k"])

rule analyze_genomes:
    input:
        aln = "work/aligned.fasta",                  # from nextalign/nextclade
        ref = "data/NC_045512.2.fasta"              # Wuhan reference
    output:
        stats = "out/stats_per_sequence.tsv",
        muts  = "out/mutations_long.tsv",
        mat   = "out/mut_matrix_sparse.npz",
        sites = "out/mut_matrix_sites.tsv",
        samp  = "out/mut_matrix_samples.tsv",
        pca   = "out/pca_kmeans.tsv",
        spike = "out/spike_aa_mutations.tsv",
        # optional (present only if k == "auto"):
        ksel  = temp("out/k_selection.tsv")
    params:
        outdir = "out",
        k      = lambda w, c, p: str(config["pca"]["k"]),
        kmin   = lambda w, c, p: int(config["pca"]["kmin"]),
        kmax   = lambda w, c, p: int(config["pca"]["kmax"])
    conda:
        "envs/analyze.yml"   # scikit-learn, biopython, numpy, pandas, scipy
    threads: 2
    shell:
        r"""
        python scripts/analyze_sarscov2.py \
            --aln {input.aln} \
            --ref {input.ref} \
            --out {params.outdir} \
            --k {params.k} \
            --kmin {params.kmin} \
            --kmax {params.kmax}

        # touch k_selection.tsv only if not created (e.g., fixed k)
        if [ ! -f {params.outdir}/k_selection.tsv ]; then
            : > {params.outdir}/k_selection.tsv
        fi
        """


