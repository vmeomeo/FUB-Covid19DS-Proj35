OUT = config["outdir"]
FIG = config["figdir"]

rule visualize:
    input:
        stats = f"{OUT}/stats_per_sequence.tsv",
        muts  = f"{OUT}/mutations_long.tsv",
        pca   = f"{OUT}/pca_kmeans.tsv",
        spike = f"{OUT}/spike_aa_mutations.tsv"
    output:
        f1 = f"{FIG}/hist_mutation_load.png",
        f2 = f"{FIG}/pca_clusters.png",
        f3 = f"{FIG}/hist_ambiguous_bases.png",
        f4 = f"{FIG}/top_spike_aas.png",
        f5 = f"{FIG}/snp_site_frequency.png",
        f6 = f"{FIG}/spike_snp_site_frequency.png"
    params:
        figdir = FIG,
        indir  = OUT
    conda: "../envs/nextalign.yaml"
    shell:
        r"""
        mkdir -p {params.figdir}
        python scripts/visualize_sarscov2.py --in {params.indir} --figdir {params.figdir}
        test -s {output.f1} || (echo 'Visualization failed' >&2; exit 1)
        """

