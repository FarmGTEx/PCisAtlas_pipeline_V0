# FASTQ Processing Pipeline

This folder contains scripts for processing raw single-cell ATAC-seq sequencing data — from demultiplexing through fragment generation and quality control. A demo dataset is included to illustrate the workflow.

## System Requirements

- **Operating System:** Linux (tested on CentOS 7/8)
- **Software Dependencies:** See `demo/environment_list.txt` for the full conda environment. Key tools include:
  - `fastp` (≥0.23) — FASTQ quality filtering
  - `chromap` (≥0.2) — read alignment for scATAC-seq
  - `samtools` (≥1.17) — BAM processing
  - `bedtools` (≥2.31) — genomic interval operations
  - `snapatac2` (≥2.3) — TSS enrichment calculation
- **Reference Files** (must be prepared separately):
  - Pig reference genome (FASTA + chromap index)
  - Gene annotation (GTF)
  - Autosome list, repeat regions, chromosome sizes, 5kb bin file

---

## Demo: Instructions to Run on Data

The `demo/` directory contains a small synthetic dataset to verify the pipeline:

| File | Description |
|------|-------------|
| `demo_I1.fastq.gz` | Index read 1 (i5 barcode) |
| `demo_I2.fastq.gz` | Index read 2 (i7 barcode) |
| `demo_R1.fastq.gz` | Read 1 (biological) |
| `demo_R2.fastq.gz` | Read 2 (biological) |
| `barcodes` | Sample-to-barcode mapping (Tissue → T5/T7/I5/I7 ranges) |
| `I5`, `I7`, `T5`, `T7` | Barcode sequence files for each index |
| `spike-in` | Spike-in barcode sequence |

### Step 1: Demultiplexing

Edit `01.demultiplex_fastq.sh` and set the demo input files, then run:

```bash
bash 01.demultiplex_fastq.sh
```

This script performs two operations:
1. `preprocess_barcodes.indel.pl` — Generates corrected barcode reference files from the raw barcode lists
2. `get_fixed.fq.indel.pl` — Demultiplexes the pooled FASTQ into per-sample paired-end FASTQ files with cell barcodes embedded in read names

**Expected output** (in `BATCH1/`):
- `BATCH1.demulti.{SampleName}_R1.fastq` and `BATCH1.demulti.{SampleName}_R2.fastq` — demultiplexed FASTQ files for each sample
- A log file recording the number of reads assigned to each sample

### Step 2: Fragment Generation

Edit `02.fastq2fragment.sh` to set the input FASTQ names, reference paths, and output prefix, then run:

```bash
bash 02.fastq2fragment.sh
```

This script performs:
1. **FASTQ filtering** (`fastp`) — adapter trimming, quality filtering (Q≥20), length filtering (≥20bp)
2. **Alignment** (`chromap`) — scATAC-seq optimized mapping to reference genome
3. **Fragment processing** (`samtools` + `perl`) — BAM → fragment file conversion, duplicate removal, mitochondrial read filtering
4. **Quality metrics** (`snapatac2`) — TSS enrichment score calculation per barcode

**Expected output**:
- `${NAME}.fragments.tsv.gz` — Fragment file (chrom, start, end, barcode, count)
- `${NAME}.metadata.txt` — Per-barcode metadata (sample, tissue, total fragments, TSS enrichment, etc.)
- `${NAME}.fastp.json` — Fastp quality report

### Step 3: Fragment QC

The `fragment.qc.{1-4}.pl` scripts apply additional cell-level filters to the fragment file:

| Script | Filter |
|--------|--------|
| `fragment.qc.1.pl` | Fragment length < 1000 bp, remove MT reads |
| `fragment.qc.2.pl` | Fragment count per barcode: > 1000 and < 100,000 |
| `fragment.qc.3.pl` | TSS enrichment > 2 |
| `fragment.qc.4.pl` | Nucleosome signal < 4, FRiP > 0.2 |

**Expected output**: `${NAME}.fragments.filtered.tsv.gz` — Filtered fragment file ready for ArchR import.

### Step 4: Visualization

```bash
Rscript QC.plot.r           # Fragment statistics and QC metrics summary
Rscript QC.plot.fraglen.r   # Fragment length distribution
Rscript Stat.fragments.r    # Per-sample fragment count summary
```

---

## Instructions for Use: How to Run on Your Data

### 1. Prepare Reference Files

Before running the pipeline, prepare the following reference files for your genome build (examples given for pig Sscrofa 11.1):

- **Genome FASTA + chromap index**: `pig.fa`, `pig.chromapindex`
- **Gene annotation**: `Sus_scrofa.Sscrofa11.1.101.chr.gtf`
- **Autosome list**: `pig_autosomes.txt` (one chromosome name per line, e.g., `chr1`, `chr2`, ...)
- **Repeat regions**: `pig.repeat.txt` (BED format)
- **5 kb bins**: `pig.5K` (BED format, 5 kb windows genome-wide)
- **Chromosome sizes**: `pig.sizes` (two-column: chr name, size in bp)

### 2. Prepare Barcode Files

Create the following files for your experiment:

- **`barcodes`**: Sample-to-index mapping. Format: `{SampleName}\t{T5_range}\t{T7_range}\t{I5_range}\t{I7_range}` (one line per sample)
- **`I5`, `I7`, `T5`, `T7`**: One barcode sequence per line for each index type
- **`spike-in`**: Spike-in barcode sequences (one per line), if applicable

The `barcodes` file uses a compact range notation: e.g., `T5-1:T5-8` means T5 lines 1 through 8 from the T5 file.

### 3. Configure and Run

1. **Edit `01.demultiplex_fastq.sh`**: Set `FQ_I1`, `FQ_I2`, `FQ_R1`, `FQ_R2` to your raw FASTQ paths, and `output` to your batch name
2. **Edit `02.fastq2fragment.sh`**: Set `NAME` to your sample/batch name, and update all `~/refgenome/` paths to your reference files. Replace `${NAME}=SAMPLE` with your sample name and iteratively run for every sample
3. Run the scripts in order:
   ```bash
   bash 01.demultiplex_fastq.sh
   bash 02.fastq2fragment.sh
   perl fragment.qc.1.pl  # ... (see individual scripts for usage)
   perl fragment.qc.2.pl
   perl fragment.qc.3.pl
   perl fragment.qc.4.pl
   ```

### 4. Verify Output

The final filtered fragment file and metadata file can be imported directly into ArchR (see `02.initialize_ArchR/`). Expected output sizes will vary by experiment scale; for reference, our atlas processed ~37 pigs across 55 tissues, yielding ~2 million nuclei after all QC filters.
