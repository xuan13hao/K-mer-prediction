# K-mer-prediction

Code and processed data for our study on predicting cancer immunotherapy
response from shotgun **metagenomic profiles** of the gut microbiome. We
compare three kinds of feature profiles built from the same sequencing reads —
**taxonomic** (MetaPhlAn), **functional / Pfam** (UProC), and raw **k-mer
counts** (KMC) — and train machine-learning classifiers to predict clinical
response (and progression-free survival, PFS) across three independent
non-small-cell lung cancer (NSCLC) cohorts.

> This repository accompanies the publication *<ADD PAPER TITLE / DOI HERE>*.

---

## Repository structure

```
K-mer-prediction/
├── README.md                 # this file
├── LICENSE                   # MIT (edit the copyright holder)
├── requirements.txt          # Python dependencies
├── preprocess/               # Step 1 — raw reads → feature tables
│   ├── README.md             #   detailed preprocessing/annotation protocol
│   ├── download.sh           #   download raw reads from SRA
│   ├── annotate.sh           #   MetaPhlAn / UProC / KMC annotation pipeline
│   ├── fq2fa.py
│   ├── removeReadSAM.py      #   remove human-mapped reads
│   ├── getMTPTable.py        #   merge MetaPhlAn output  → otu.count.csv
│   ├── getUproCTable.py      #   merge UProC output      → pfam.count.tsv
│   ├── getKMCTable.py        #   merge per-sample k-mer counts
│   ├── filterKMCTable.py     #   prevalence-filter the k-mer table
│   └── mergeKMCTable/        #   C++ external-merge tool (main.cpp + makefile)
├── analysis/                 # Step 2 — machine-learning prediction
│   ├── predict.py            #   train MLP, output per-iteration probabilities
│   ├── reproduce.py          #   re-run a single fixed-seed configuration
│   ├── compute_importance.py #   Random-Forest feature importance
│   └── run_experiments.sh    #   driver for the main paper experiments
├── figures/                  # Step 3 — figures & tables
│   └── make_figures.ipynb    #   reads data/results/*, draws the bar plots
└── data/
    ├── README.md             # data availability + expected file layout
    └── results/              # (empty) prediction outputs are written here
```

**The processed datasets are NOT stored in this repository** (to keep it
lightweight). Download them from the link below and unzip into `data/`:

> **https://cbb.ittc.ku.edu/software/kmer/data.zip**

See [data/README.md](data/README.md) for the file list and layout.

---

## Installation

The analysis was run with **Python 3.7.9**.

```bash
pip install -r requirements.txt
```

To reproduce the feature tables from raw sequencing reads you additionally need
the following external bioinformatics tools (see
[preprocess/README.md](preprocess/README.md) for versions and exact commands):

FastQC, Trimmomatic, Fiona, BWA, MetaPhlAn, UProC, KMC3, and (optionally) the
SRA Toolkit.

---

## Data

The processed datasets (sample metadata, taxonomic / Pfam / KEGG / k-mer count
matrices, and pre-computed prediction outputs) are **not stored in git** because
of their size. Download the archive and unzip it into `data/`:

> **https://cbb.ittc.ku.edu/software/kmer/data.zip**

(Paste the link into any browser.) Raw sequencing reads are also available from
NCBI (BioProjects below). See [data/README.md](data/README.md) for the full file
list and the expected layout.

---

## Reproducing the results

### Step 1 — Preprocessing (optional)

Only needed if you want to rebuild the feature tables from raw reads. Raw data
are available from NCBI:

| Cohort | BioProject  |
|--------|-------------|
| DS1    | PRJNA866654 |
| DS2    | PRJEB22863  |
| DS3    | PRJNA751792 |

Follow [preprocess/README.md](preprocess/README.md) for download, quality
control, human-read removal, error correction, annotation (MetaPhlAn / UProC /
KMC), and table merging.

### Step 2 — Prediction

`analysis/predict.py` trains an MLP classifier on a chosen feature profile. For
each of 100 iterations it searches up to 1000 random seeds for the best model
(selected by F1 on a 70/30 stratified split) and appends the test-set
probabilities to `data/results/<profile>_<dataset>_<report_name>_prob.csv`.

```bash
cd analysis
# <profile> <report_name> <sp_code> <feature1> <feature2> ...
python3 predict.py pfam methanogenesis mix PF01993 PF02240 PF02745   # example
```

- `<profile>`: `pfam`, `tax`, `k16`–`k25`, `K30`
- `<sp_code>`: `mix`, `dds3`, `ds3`, `ds12`, `ds12t3`, `ds3t12`, `ds123t4`,
  `ds4t123`, `ds123t5`, `ds5t123` (train/test split — run `python3 predict.py`
  with no arguments for the full description)
- features: specific feature IDs, or `all`

Run all the main experiments at once:

```bash
cd analysis
bash run_experiments.sh
```

Feature importance (Random Forest, 100 iterations) →
`data/imp/RF_<profile>_<sp_code>_importance.csv`:

```bash
python3 compute_importance.py pfam mix
```

Verify a single configuration with fixed seeds (prints confusion matrix, F1,
macro-F1, MCC):

```bash
python3 reproduce.py <rs_model> <rs_sp> <sp_code> <profile> <feature1> ...
```

### Step 3 — Figures

Open the notebook and run all cells. It reads the probability files in
`data/results/` and reproduces the paper's bar plots and P-value tables. The
result files are produced by Step 2 (`run_experiments.sh`); they are not shipped
with the repository.

```bash
cd figures
jupyter notebook make_figures.ipynb
```

Significance markers on the plots: `*` P<5e-3, `**` P<5e-4, `***` P<1e-5-ish
bracket, `****` P<1e-5 (see the thresholds in the notebook).

---

## Script reference

| Script | Purpose |
|--------|---------|
| `preprocess/download.sh` | Download raw FASTQ from SRA |
| `preprocess/annotate.sh` | Per-sample MetaPhlAn / UProC / KMC annotation |
| `preprocess/getMTPTable.py` | Merge MetaPhlAn profiles → `otu.count.csv` |
| `preprocess/getUproCTable.py` | Merge UProC profiles → `pfam.count.tsv` |
| `preprocess/getKMCTable.py` / `mergeKMCTable/` | Merge per-sample k-mer counts |
| `preprocess/filterKMCTable.py` | Prevalence-filter the k-mer table |
| `analysis/predict.py` | Train MLP, output prediction probabilities |
| `analysis/compute_importance.py` | Random-Forest feature importance |
| `analysis/reproduce.py` | Re-run one fixed-seed configuration |
| `analysis/run_experiments.sh` | Driver for the main experiments |
| `figures/make_figures.ipynb` | Generate figures and P-value tables |

---

<!-- ## Citation

If you use this code or data, please cite:

```
<ADD BIBTEX / CITATION HERE>
``` -->

## License

Released under the MIT License — see [LICENSE](LICENSE).
