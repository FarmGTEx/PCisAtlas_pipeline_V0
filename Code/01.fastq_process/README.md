🔬# Demo Sequencing Data for Demultiplexing
This folder includes a demo dataset and associated scripts to illustrate the demultiplexing and preprocessing workflow.

💻## System Requirements
Operating System:
Linux (tested on: Linux login01 3.10.0-862.el7.x86_64)

Software Dependencies:
Please refer to the file environment_list.txt in this directory for a full list of required software and version information.

▶️## Running the Demo
Step 1: Demultiplexing
Run the script:
bash 01.demultiplex_fastq.sh

This step will demultiplex the raw sequencing data and generate paired FASTQ files for each sample. Corrected cell barcodes will be added to the read names.

Step 2: Preprocessing
Run the script:
bash 02.fastq2fragment.sh

This step will perform:
FASTQ filtering
Alignment
Fragment filtering and reformatting

The final output includes:
A fragment file
Related metadata files
