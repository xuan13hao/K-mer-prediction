# `preprocess/` — from raw reads to feature tables

This folder documents how the raw metagenomic sequencing reads were turned into
the taxonomic / Pfam / k-mer count tables used for prediction. You only need
this if you want to rebuild the tables from scratch; otherwise download the
processed tables (see [../data/README.md](../data/README.md)).

## 1. Download

If you have the [SRA Toolkit](https://github.com/ncbi/sra-tools) installed, use
`download.sh` (edit the `fastq-dump` path at the top to point at your
`sratoolkit/bin/`). Otherwise download the raw data however you prefer using the
BioProject accessions:

| Cohort | BioProject  | Link |
|--------|-------------|------|
| DS1 | PRJNA866654 | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA866654 |
| DS2 | PRJEB22863  | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJEB22863  |
| DS3 | PRJNA751792 | https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA751792 |

Notes:
- Not all DS2 samples are NSCLC-related; we only used the samples taken at
  timepoint **T0**.
- Hosting the raw data requires **> 2 TB** of disk.

## 2. Quality control, host removal, error correction

Required tools (install and set each tool's path):
FastQC (v0.11.8), Trimmomatic (0.38), Fiona (0.2.10), BWA (0.7.16a),
MetaPhlAn (4.0.3), UProC (1.2.0), KMC3 (3.2.1).

**Quality control** — inspect with FastQC, then remove adapters / low-quality
bases with Trimmomatic:

```bash
# DS1 (paired-end, known adapters)
java -jar trimmomatic-0.38.jar PE $ID.read1.fastq $ID.read2.fastq \
    $ID.read1.paired.fq $ID.read1.unpaired.fq $ID.read2.fq $ID.read2.unpaired.fq \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:35 -threads 16

# DS2 / DS3 (single-end)
java -jar trimmomatic-0.38.jar SE $ID.fastq $ID.trimed.fastq \
    ILLUMINACLIP:TruSeq2-SE.fa:2:30:10 CROP:280 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:35 -threads 16
```

**Remove human reads** — map to hg38 with BWA and drop mapped reads with
`removeReadSAM.py`:

```bash
bwa index hg38.fasta
bwa mem -t 16 hg38.fa $ID.trimed.fastq > $ID.sam
python3 removeReadSAM.py $ID.trimed.fastq $ID.sam > $ID.trimed.clean.fq
```

**Error correction** (recommended):

```bash
fiona -nt 16 -g 4639675 $ID.trimed.clean.fq $ID.trimed.clean.ec.fq
```

## 3. Annotation

```bash
bash annotate.sh "$FILENAMES" $SAMPLE_ID
```

- `$FILENAMES` is `"$ID.read1.paired.clean.fq $ID.read2.paired.clean.fq"` for
  DS1 (keep the quotes!) and `$ID.trimed.clean.fq` for DS2/DS3.
- `$SAMPLE_ID` is the label for this sample (e.g. its SRA accession).
- Produces three folders per sample: `otu/` (MetaPhlAn taxonomy),
  `pfam/` (UProC functional), and `K30/` (30-mer counts).
- The 30-mer counts of all 417 samples require **> 2 TB** of disk.

## 4. Merge into OTU-like tables

```bash
# Taxonomic table  → otu.count.csv
python3 getMTPTable.py metadata.csv otu

# Functional (Pfam) table → pfam.count.tsv
python3 getUproCTable.py metadata.csv pfam
```

Both scripts need `pandas`. Run either with no arguments to see its help.

**k-mer table** — merge the per-sample KMC outputs with the C++ tool in
`mergeKMCTable/`:

```bash
cd mergeKMCTable && make          # needs the boost library
# move the built binary + KMC_list.txt into the K30/ folder, then:
./getKMCTable KMC_list.txt        # → k30.count.tsv
```

`mergeKMCTable` implements a two-way external merge; tune it to your available
memory/disk. `getKMCTable.py` is a Python helper for the same merge step.

**Prevalence-filter the k-mer table** (optional) with `filterKMCTable.py`:

```bash
# <metadata.csv> <inName> <outName> [<ratio percent>]
python3 filterKMCTable.py ../data/metadata.csv k30.count.tsv k30.filtered 60
```

The metadata path is now a command-line argument (previously hard-coded).

---

See the paper's supplementary material for the full rationale behind each step.
