DATASET = config["dataset_fasta"]
REF     = config["reference_fasta"]
WORK    = config["workdir"]

rule align_nextalign:
    input:
        fasta = DATASET,
        ref   = REF
    output:
        aligned = f"{WORK}/aligned.fasta"
    conda: "../envs/nextalign.yaml"   # note: path is relative to this rules/ file
    threads: 2
    shell:
        r"""
        mkdir -p {WORK}
        # New nextalign CLI: positional input FASTA, -r for reference, -O for output dir
        nextalign run -r {input.ref} -O {WORK} {input.fasta}

        # nextalign writes something like: {WORK}/nextalign.aligned.fasta (or similar).
        # Find it and normalize name to aligned.fasta for downstream rules.
        ALN=$(ls {WORK}/*aligned*.fasta 2>/dev/null | head -n1)
        test -s "$ALN" || (echo "No aligned FASTA produced by nextalign" >&2; exit 1)
        cp "$ALN" {output.aligned}
        """

