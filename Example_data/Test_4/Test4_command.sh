#!/bin/bash

which GALEON_ControlScript.py 

# Run Galeon
GALEON_ControlScript.py clusterfinder -a GFFs/ -g 100,300 -e enabled -p Proteins/ -pm True -outdir 4_OneFam_PhysEvoDistances_GFF3_pm/

# Get evolutionary statistics (Cst) and perform the Mann-Whitney test
GALEON_GetEvoStats.py -clust 4_OneFam_PhysEvoDistances_GFF3_pm/ -prot Proteins/ -coords GFFs

# Generate summary files and tables for GR family
GALEON_SummaryFiles.py -fam GR -clust 4_OneFam_PhysEvoDistances_GFF3_pm/ -coords GFFs -ssize ChrSizes.txt -sfilter 7

# Generate summary files and tables for IR family
GALEON_SummaryFiles.py -fam IR -clust 4_OneFam_PhysEvoDistances_GFF3_pm/ -coords GFFs -ssize ChrSizes.txt -sfilter 7

# Create a summary report
GALEON_Report.py -clust 4_OneFam_PhysEvoDistances_GFF3_pm/ -plots Plots -ssize ChrSizes.txt -echo False