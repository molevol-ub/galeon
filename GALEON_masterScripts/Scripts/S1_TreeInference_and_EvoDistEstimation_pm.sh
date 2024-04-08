#!/bin/bash


# | 1 | # Inputs
Prot_dir=$1 # Input directory with proteins to be processed by this pipeline (*.fasta and *.fa extensions are admitted)

# | 2 | # Parameteres
software=$2 # Program to infer the phylogenetic tree across the input proteins (FastTree or iqtree)

# master script 
# run_evodistance_fromaln.sh to be run  (don't modify it)
evodist_pipe=$3 

# Bin directory
bin_dir=$4


# ~~~~ Execution ~~~~~ #

# Count how many FASTA files are there
# Use ls and grep to count files with .txt extension
count=$(ls "$Prot_dir"/*.aln 2>/dev/null | grep -c .)

# Check if count is zero
if [ "$count" -eq 0 ]; then
  echo "Error: No MSA files found in '$Prot_dir' input directory."
  exit 1  # Exit the script with an error code
else
  echo "- Number of MSA files in '$Prot_dir' input directory: $count"
fi

# Run the master script "run_evodistance.sh" over each fasta file
# - using "mafft-linsi" aligment algorithm 
# - and either "FastTree" or "iqtree" as program to infer the phylogenetic tree across the input proteins

for file in $Prot_dir/*aln
do
	f_extension="${file##*.}"
	f_basename="$(basename -- $file)"
	f_name="${f_basename%.*}"

    echo $file
    echo ${Prot_dir}/$f_basename
    bash $evodist_pipe -pm ${Prot_dir}/$f_basename -t $software -o ${Prot_dir}/$f_name -b $bin_dir

done
