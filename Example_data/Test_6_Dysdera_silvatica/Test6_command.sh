#!/bin/bash

which GALEON_ControlScript.py 

# Run Galeon
GALEON_ControlScript.py clusterfinder -a GFFs/ -g 100 -e enabled -p Proteins/ -outdir 6_OneFam_PhysEvoDistances_GFF3/ -f orange -pm True

# Get evolutionary statistics (Cst) and perform the Mann-Whitney test
GALEON_GetEvoStats.py -clust 6_OneFam_PhysEvoDistances_GFF3/ -prot Proteins/ -coords GFFs

# Generate summary files and tables for GR family
GALEON_SummaryFiles.py -fam GR -clust 6_OneFam_PhysEvoDistances_GFF3/ -coords GFFs -ssize ChrSizes.txt -sfilter 7 -colval magma

# Generate summary files and tables for IR family
GALEON_SummaryFiles.py -fam IR -clust 6_OneFam_PhysEvoDistances_GFF3/ -coords GFFs -ssize ChrSizes.txt -sfilter 7 -colval magma

# Create a summary report
GALEON_Report.py -clust 6_OneFam_PhysEvoDistances_GFF3/ -plots Plots -ssize ChrSizes.txt -echo False
