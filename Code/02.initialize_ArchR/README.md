# ArchR Analysis
This folder contains scripts to initialize an ArchR project and perform downstream analyses, including matrix normalization, dimensionality reduction, clustering, and peak calling.

## System Requirements
Operating System:
Linux (tested on: Linux 4.18.0-513.5.1.el8_9.x86_64)

## Software Dependencies:
Please refer to the files environment_list.txt and r_packages_version.txt in this directory for the complete list of required software and R package versions.

## Included Scripts
01.create_ArchR_object.r: Initializes the ArchR project and imports input files.
02.initially_clustering.r: Performs quality control, dimensionality reduction (e.g., LSI), and unsupervised clustering.
03.call_peaks.r: Conducts peak calling based on clustered data.
