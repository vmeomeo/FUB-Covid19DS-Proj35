REF = config["reference_fasta"]

rule get_reference:
    output: REF
    conda: "../envs/nextalign.yaml"
    shell:
        r"""
        mkdir -p $(dirname {output})
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text" > {output}
        test -s {output} || (echo 'Reference download failed' >&2; exit 1)
        """

