## **Project 3: Tasks**
**Project 3: Development of a risk score to analyze SARS-CoV-2 genomes**

- Download SARS-CoV-2 sequence set from https://osf.io/vd7gj
- Align to the Index SARS-CoV-2 sequence from the Wuhan outbreak (NCBI: NC_045512.2)
- Find mutant positions (e.g. CIGAR string or via some tool)
- Find mutant positions in the Spike gene
- (can do in sample in project 5) Perform basic statistics
- Number of ambiguous nucleotides per sequence
- Number of mutations (SNPs, INDELs)
- Project and cluster sequences (for example: PCA, kmeans)
- (optional) Download SARS-CoV-2 antibody escape scores provided by the  Bloom Lab (paper: https://www.science.org/doi/10.1126/science.abf9302, online resource: https://github.com/jbloomlab/SARS2_RBD_Ab_escape_maps)
- Develop a "risk" score, e.g. based on the Spike sequences and the Bloom scores (Can skip)
- Among your example data, what's the most dangerous SARS-CoV-2 sample in your eyes? (Can skip)

### Directory map
```
project3/
├── data/
│   ├── sequences.fasta              # put your OSF fasta here
│   └── ref/NC_045512.2.fasta        # auto-downloaded by rule
├── envs/
│   ├── base.yaml
│   └── minimap2.yaml
├── workflow/
│   ├── Snakefile
│   └── scripts/
│       ├── parse_minimap2_cs.py
│       └── pca_kmeans_and_risk.py
└── out/                             # results land here
```
You can `mkdir <directory_name>` to create the mentioned directory if you haven't had it.

### Step 1: Setup
Check `setup.MD` at the parent directory

### Step 2: Data
Run in terminal (Debian/Ubuntu):

`cd project3/`

`wget --no-check-certificate https://osf.io/vd7gj/download -O data/risk-assessment-sequences.fasta`

Get reference FASTA via NCBI E-utilities (plain FASTA):

```mkdir ref
wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text" -O ref/NC_045512.2.fasta```



