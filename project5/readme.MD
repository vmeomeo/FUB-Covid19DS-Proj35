# **Project 5: Cross-Immunization**

- Get a GISAID user account: [GISAID] https://gisaid.org/register/ This will take a few days! → Done
- Write a small R script to download mutation profiles from outbreak.info API for the following SARS-CoV-2 variants using your gisaid account:
    - JN.1
    - JN.2
    - JN.3
    - KP.3
    - XBB.1.5
- Use the methods described in the Hands-on to obtain the mutation profiles of these variants and compute their cross-reactivity. Plot the result in a cross-reactivity map.

- Reproduce the Relative fitness Trends of 5 of the variants (of your choice) shown in the Figure 3 of Manuscript below:
    - Manuscript: [SARS-CoV2 Evolution on a dynamic Immune landscape](https://www.nature.com/articles/s41586-024-08477-8)
    - It does not need to look exactly the same, but you should be able to reproduce the trends.
    - NB: the appendage .X refers to all daughter variants with the parent included, and the trend is re-normalized to their respective proportions (see Hints: Project5_hints.ipynb found [here](https://github.com/AlexiaNomena/SC2_VASIL/))

- Analyze results.
- Prepare a talk to present the results & analysis.

## Directory map
```
project3/
├── data/
│   ├── dms_per_ab_per_site.csv      # provided
│   ├── PK_for_all_Epitopes.csv      # From other publications
│   ├── mut_profiles                 # GISAID mutation profiles (will be downloaded by user)
│   └── input/test_Germany           # Germany genomic surveillance data (provided)
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
├── logs/                            # snakemake log files
│
├── figures/                         # plots generate in the notebook will be here
│
└── out/                             # results land here
```
You can `mkdir <directory_name>` to create the mentioned directory if you haven't had it.

## How to use this repo?
To download the data, make sure your current working directory is `FUB-Covid19DS-Proj35/project5`.

How to check that? See the path in the terminal, or type `pwd` command. The result should be `{something}/FUB-Covid19DS-Proj35/project5`

Make sure you activate the mamba environment before running (`mamba activate envs/workshop`)
Then you run `snakemake --cores 4` to download mutation profiles. You can edit the list of lineages at the `config/config.yaml`

After that, you can run the jupyter notebook `workflow/notebook/notebook.ipynb` to analyse the data and create plots.

# Citation
Raharinirina, N. A., Gubela, N., Börnigen, D., Smith, M. R., Oh, D. Y., Budt, M., ... & von Kleist, M. (2025). SARS-CoV-2 evolution on a dynamic immune landscape. Nature, 639(8053), 196-204.
