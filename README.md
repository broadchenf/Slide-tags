# Slide-tags
This guide explains how to run our spatial positioning pipeline assigning coordinates to profiled nuclei. We assume you have (1) demultiplexed your spatial barcode library sequencing data into R1 and R2 fastq files, and (2) run your gene expression data through Cell Ranger (and CellBender if necessary). We are working on a more user-friendly pipeline that will be available in the coming months. 

#### 1. FASTQ processing (sb_processing.sh)

This script greps for the UP site in spatial barcode reads and downsamples resulting cleaned library to 25 million reads. 


#### 2. Match spatial barcode reads to 10x cell barcodes (cellbarcode_matching.R)

This script matches real 10x cell barcodes (whitelist from Cell Ranger or CellBender cell calls) with spatial barcode sequences from the spatial barcode library. 


#### 3. Assign coordinates to spatial barcodes (bead_matching.py)
This script matches candidate spatial barcodes with whitelisted cell barcodes to in situ sequenced spatial barcode coordinates. 


#### 4. Spatially map nuclei (spatial_positioning.R)
Run this script to assign spatial coordinates to profiled nuclei using DBSCAN. Several DBSCAN minPts parameters are tested under a constant eps. The parameter set with the highest proportion of positioned cells (single DBSCAN cluster) is chosen as the optimal set. 
