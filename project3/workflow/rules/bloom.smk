OUT = config["outdir"]

rule bloom_risk:
    input:
        stats = f"{OUT}/stats_per_sequence.tsv",
        muts  = f"{OUT}/mutations_long.tsv"
    output:
        risk = f"{OUT}/stats_with_bloom.tsv"
    conda: "../envs/nextalign.yaml"
    shell:
        r"""
        python scripts/compute_bloom_risk.py
        test -s {output.risk} || (echo 'Bloom risk step failed' >&2; exit 1)
        """

