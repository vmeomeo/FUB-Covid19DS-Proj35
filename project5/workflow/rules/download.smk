rule download_gisaid:
    output:
        "data/mut_profiles/{lineage}.txt"
    params:
        lineage=lambda wildcards: wildcards.lineage.replace('_', '.')
    log:
        "logs/download_gisaid_{lineage}.log"
    shell:
        """
        Rscript workflow/scripts/download_gisaid.R {params.lineage} {output} > {log} 2>&1
        """