# Snakefile or rules/analyze.smk

REF   = config["reference_fasta"]
WORK  = config["workdir"]
OUT   = config["outdir"]

rule analyze_genomes:
    input:
        aln = f"{WORK}/aligned.fasta",        # nextalign output
        ref = REF                             # Wuhan reference
    output:
        stats   = f"{OUT}/stats_per_sequence.tsv",
        muts    = f"{OUT}/mutations_long.tsv",
        mat     = f"{OUT}/mut_matrix_sparse.npz",
        sites   = f"{OUT}/mut_matrix_sites.tsv",
        samp    = f"{OUT}/mut_matrix_samples.tsv",
        pca     = f"{OUT}/pca_kmeans.tsv",
        spike   = f"{OUT}/spike_aa_mutations.tsv",
        metrics = temp(f"{OUT}/kmeans_metrics.tsv")  # optional; script creates it in auto mode
    params:
        outdir   = OUT,
        k        = lambda wc: str(config["pca"]["k"]),          # e.g. "auto" or "3"
        kmin     = lambda wc: int(config["pca"]["kmin"]),       # e.g. 2
        kmax     = lambda wc: int(config["pca"]["kmax"]),       # e.g. 10
        kselect  = lambda wc: str(config["pca"].get("select", "silhouette"))  # "silhouette", "ch", "db"
    conda:
        "../envs/nextalign.yaml"
    threads: 2
    shell:
        r"""
        python scripts/analyze_sarscov2.py \
            --aln {input.aln} \
            --ref {input.ref} \
            --out {params.outdir} \
            --k {params.k} \
            --kmin {params.kmin} \
            --kmax {params.kmax} \
            --k_select {params.kselect}

        # Ensure the optional metrics file exists so Snakemake is happy
        if [ ! -f {params.outdir}/kmeans_metrics.tsv ]; then
            : > {params.outdir}/kmeans_metrics.tsv
        fi
        """


