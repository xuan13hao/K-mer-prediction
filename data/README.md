# `data/` — data availability

The processed datasets used in this study (sample metadata, taxonomic / Pfam /
KEGG / k-mer count matrices, and the pre-computed prediction outputs) are a
**controlled research asset and are NOT distributed with this repository.**

## How to obtain the data

The data are **available from the authors upon reasonable request**. Please
contact the corresponding author (see the paper) with a short description of
your intended use. Access is granted on a per-application basis.

Raw sequencing reads can also be obtained directly from NCBI under the
following BioProjects:

| Cohort | BioProject  |
|--------|-------------|
| DS1    | PRJNA866654 |
| DS2    | PRJEB22863  |
| DS3    | PRJNA751792 |

and re-processed into feature tables using the scripts in
[../preprocess/](../preprocess/).

## Where to put the files once you receive them

After you obtain the processed tables, place them in this `data/` folder so the
scripts can find them (every script resolves this folder as `../data`). The
expected files are:

| File | Format | Description |
|------|--------|-------------|
| `metadata.csv` | CSV | Per-sample metadata: Run id, SampleID, age, BMI, sex, best response, OS/PFS, antibiotics, PD-L1, timepoint, PFS buckets, cohort. |
| `otu.count.csv` | CSV | Species-level taxonomic count matrix (rows = taxa, columns = samples). |
| `pfam.count.tsv` | TSV | Pfam (functional) raw-count matrix. |
| `kegg.count.tsv` | TSV | KEGG pathway raw-count matrix. |
| `k30.count.tsv` | TSV | 30-mer count matrix. |
| `k16.count.pvalue.tsv`, `k19_24.count_33_melanowa_PRJNA615114.csv`, ... | TSV/CSV | k-mer count tables for the corresponding `predict.py` profile options. |
| `pfam.acclengname.csv` | CSV | Pfam accession → length → name map (RPKM normalisation). |

## Table conventions

- **CSV** = comma-separated, **TSV** = tab-separated.
- Rows are features (taxa / Pfam families / k-mer ids); columns are samples;
  the first column is the feature id/index.

## Generated sub-folders (not tracked)

Running the analysis creates these (git-ignored) folders:

- `results/` — `*_prob.csv` prediction outputs from `analysis/predict.py`,
  read by `figures/make_figures.ipynb`
- `imp/` — feature-importance CSVs from `analysis/compute_importance.py`
- `reproduce/` — scratch folder created by `analysis/reproduce.py`
