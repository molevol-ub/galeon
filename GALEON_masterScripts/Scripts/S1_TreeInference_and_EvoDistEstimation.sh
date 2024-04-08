#!/bin/bash


# | 1 | # Inputs
Prot_dir=$1 # Input directory with proteins to be processed by this pipeline (*.fasta and *.fa extensions are admitted)

# | 2 | # Parameteres
software=$2 # Program to infer the phylogenetic tree across the input proteins (FastTree or iqtree)

# master script 
# run_evodistance.sh to be run  (don't modify it)
evodist_pipe=$3 

# Bin directory
bin_dir=$4


## Mafft aligment option
mafft_default="mafft-linsi"
mafft_algorithm=$5 # mafft-linsi or mafft

# Check if mafft_algorithm is provided
if [ -z "$mafft_algorithm" ]; then
    # echo "# MSA algorithm. Using default value: $mafft_default"
    mafft_algorithm=$mafft_default
fi


# ~~~~ Execution ~~~~~ #

# Count how many FASTA files there are?
# Use ls and grep to count files with .txt extension
count=$(ls "$Prot_dir"/*.fasta 2>/dev/null | grep -c .)

# Check if count is zero
if [ "$count" -eq 0 ]; then
  echo "Error: No FASTA files found in '$Prot_dir' input directory."
  exit 1  # Exit the script with an error code
else
  echo "- Number of FASTA files in '$Prot_dir' input directory: $count"
fi

# Run the master script "run_evodistance.sh" over each fasta file
# - using "mafft-linsi" aligment algorithm 
# - and either "FastTree" or "iqtree" as program to infer the phylogenetic tree across the input proteins

for file in $Prot_dir/*
do
	f_extension="${file##*.}"
	f_basename="$(basename -- $file)"
	f_name="${f_basename%.*}"


	if [[ $f_extension == "fasta" ]]
	then
		out_aln=$f_name.aln

		# echo "# ~ Input protein FASTA file: '${f_basename}'"
		# echo "# ~ Output alignment file: '${out_aln}'"

		# echo "${f_name}	${Prot_dir}/${f_name}_distancematrix.tsv" >> ../$out_evofile_list
		bash $evodist_pipe -p ${Prot_dir}/$f_basename -t $software -o ${Prot_dir}/$f_name -b $bin_dir -m $mafft_algorithm

	fi

done
	# elif [[ $f_extension == "fa" ]]
	# then
	# 	out_aln=$f_name.aln

	# 	# echo "# ~ Input protein FASTA file: '${f_basename}'"
	# 	# echo "# ~ Output alignment file: '${out_aln}'"

	# 	# echo "${f_name}	${Prot_dir}/${f_name}_distancematrix.tsv" >> ../$out_evofile_list
	# 	bash $evodist_pipe -p ${Prot_dir}/$f_basename -t $software -o ${Prot_dir}/$f_name -b $bin_dir

	# else
		# echo "# Warning! This file '${f_basename}' does not have a known FASTA extension (*.fasta, *.fa). Nothing is done here..."

# 	fi

# done
