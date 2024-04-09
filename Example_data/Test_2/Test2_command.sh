#!/bin/bash

which GALEON_ControlScript.py 

# Run Galeon
GALEON_ControlScript.py clusterfinder -a GFFs/ -g 100,200 -e disabled -outdir 2_OneFam_PhysDistOnly_GFF3

# Generate summary files and tables for GR family
GALEON_SummaryFiles.py -fam GR -clust 2_OneFam_PhysDistOnly_GFF3/ -coords GFFs -ssize ChrSizes.txt -sfilter 7

# Generate summary files and tables for IR family
GALEON_SummaryFiles.py -fam IR -clust 2_OneFam_PhysDistOnly_GFF3/ -coords GFFs -ssize ChrSizes.txt -sfilter 7

# Create a summary report
GALEON_Report.py -clust 2_OneFam_PhysDistOnly_GFF3/ -ssize ChrSizes.txt -echo False