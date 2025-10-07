REF   = config["reference_fasta"]
WORK  = config["workdir"]
OUT   = config["outdir"]
K     = int(config["kmeans_k"])

rule analyze:
    input:
        aligned = f"{WORK}/aligned.fasta",
        ref     = REF
    output:
        stats   = f"{OUT}/stats_per_sequence.tsv",
        muts    = f"{OUT}/mutations_long.tsv",
        pca     = f"{OUT}/pca_kmeans.tsv",
        spikeaa = f"{OUT}/spike_aa_mutations.tsv",
        matrix  = f"{OUT}/mut_matrix_sparse.npz"
    params:
        outdir = OUT,
        k = K
    conda: "../envs/nextalign.yaml"
    shell:
        r"""
        mkdir -p {params.outdir}
        python scripts/analyze_sarscov2.py \
          --aln {input.aligned} \
          --ref {input.ref} \
          --out {params.outdir} \
          --k {params.k}
        test -s {output.stats} || (echo 'Analysis failed' >&2; exit 1)
        """

