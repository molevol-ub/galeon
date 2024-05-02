#!/bin/bash

which GALEON_ControlScript.py 

# Run Galeon
GALEON_ControlScript.py clusterfinder -a GFFs/ -F BetweenFamilies -e disabled -outdir 5_TwoFamJointAnalysis_GFF3

# Generate summary files and tables
GALEON_SummaryFiles.py -fam merged -clust 5_TwoFamJointAnalysis_GFF3/ -coords GFFs/merged_dir/ -ssize ChrSizes.txt -sfilter 7

# Create a summary report
GALEON_Report.py -clust 5_TwoFamJointAnalysis_GFF3/ -ssize ChrSizes.txt -echo False
