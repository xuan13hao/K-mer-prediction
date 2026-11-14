#!/usr/bin/env bash
#
# run_experiments.sh
# -------------------
# Reproduce the main prediction experiments reported in the paper.
#
# Every command below calls `predict.py`, which trains an MLP classifier on the
# requested feature profile and writes its per-iteration prediction
# probabilities to:
#
#     ../data/results/<profile>_<dataset>_<report_name>_prob.csv
#
# Those *_prob.csv files are exactly what `../figures/make_figures.ipynb`
# reads to draw the figures and compute the P-values.
#
# The processed datasets required by predict.py are NOT shipped with this
# repository (see ../data/README.md). Obtain them and place them in ../data/
# before running this script.
#
# Usage:
#     cd analysis
#     bash run_experiments.sh
#
# NOTE: each command performs a 1000-seed random search over 100 iterations and
# is therefore CPU-intensive (hours per command). Run selectively if needed.

set -euo pipefail
cd "$(dirname "$0")"

# Datasets / train-test splits (see predict.py help for the full list):
#   mix    [1,2,3 | 1,2,3]   all samples from all three datasets
#   dds3   [1,2,3 | 3]       train on DS1+DS2 plus part of DS3
#   ds3    [3 | 3]
#   ds12   [1,2 | 1,2]
DATASETS=(mix dds3 ds3 ds12)

# ---------------------------------------------------------------------------
# 1. Functional (Pfam) feature sets: one-carbon metabolism & methanogenesis
# ---------------------------------------------------------------------------
ONE_CARBON="PF05221 PF02741 PF01913 PF02007"
METHANOGENESIS="PF01993 PF02240 PF02745 PF02662 PF02505 PF03201 PF04211 PF04210 PF05440 PF02249 PF06253 PF09505 PF09472"
ONE_CARBON_METHANOGENESIS="$ONE_CARBON $METHANOGENESIS"

for ds in "${DATASETS[@]}"; do
    python3 predict.py pfam one_carbon              "$ds" $ONE_CARBON
    python3 predict.py pfam methanogenesis          "$ds" $METHANOGENESIS
    python3 predict.py pfam 1C_methanogenesis       "$ds" $ONE_CARBON_METHANOGENESIS
done

# ---------------------------------------------------------------------------
# 2. Taxonomic feature: Akkermansia muciniphila
# ---------------------------------------------------------------------------
for ds in "${DATASETS[@]}"; do
    python3 predict.py tax Akkermansia_muciniphila  "$ds" Akkermansia_muciniphila
done

# ---------------------------------------------------------------------------
# 3. Cross-dataset generalisation (train on one cohort, test on another)
#    ds12t3 = train on DS1+DS2, test on DS3
#    ds3t12 = train on DS3,      test on DS1+DS2
# ---------------------------------------------------------------------------
python3 predict.py pfam ds12t3_pfam ds12t3 all
python3 predict.py tax  ds12t3_tax  ds12t3 all
python3 predict.py pfam ds3t12_pfam ds3t12 all
python3 predict.py tax  ds3t12_tax  ds3t12 all

echo "All experiments finished. Results written to ../data/results/"
