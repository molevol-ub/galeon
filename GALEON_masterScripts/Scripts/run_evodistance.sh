#!/bin/bash

### Script to retrieve the evolutionary distance between proteins

## Requires the following binaries (included in bin)
## Mafft
# Fasttree
# Newick utils: https://bioweb.pasteur.fr/packages/pack@newick-utils@1.6

VERSION=1.0

SCRIPT=$(realpath "$0")
SCRIPTPATH=$(dirname "$SCRIPT")

PROTFILE=''
ALNMETHOD=mafft-linsi
TREEMETHOD=FastTree
BINFOLDER="$SCRIPTPATH/bin"
OUTPUT=Output

function usage {
	echo
	echo "Usage:"
	echo "  $(basename $0)"
	echo "      -p Input protein fasta file"
	echo "      -m Program to conduct the multiple sequence alignment. Specify 'mafft-linsi' or 'mafft', mafft-linsi is more accurate but slow for hundreds of protein sequences (Default=$ALNMETHOD)"	
	echo "      -t Program to infer the phylogenetic tree across the input proteins. Specify 'FastTree' or 'iqtree' (Default=$TREEMETHOD)"
	echo "      -b PATH to bin folder containing the required programs (Default=$BINFOLDER)"
	echo "      -o Prefix for the generated output files (Default=$OUTPUT)"
	echo
	exit 0
}

# Read options

if [ "$#" -lt "1" ]; # at least 1 argument
then
	usage
fi

while [ $# -gt 0 ]; do
	case "$1" in
		-h|-help) usage
				;;
		-p) shift
			PROTFILE=$1
			;;				
		-m) shift
			ALNMETHOD=$1
			;;
		-t) shift
			TREEMETHOD=$1
			;;
		-b) shift
			BINFOLDER=$1
			;;		
		-o) shift						# *mod*
			OUTPUT=$1					# *mod*
			;;							# *mod*
		*)	echo 
			echo "ERROR - Invalid option: $1"
			echo
			usage
			;;
	esac
	shift
done

# Start script

export PATH=$BINFOLDER:$PATH


# Run alignment
# echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "# | S1.1 Step |"
echo -e "# - Multiple Sequence Alignment - #\n"
echo -e "Start\n"


if [ $ALNMETHOD == "mafft-linsi" ] ; then

	if ! command -v mafft-linsi &> /dev/null
	then	
		echo -e "Cannot find mafft-linsi binary, make sure to export the path";
		usage
		exit 1;
	fi

	echo " - Input protein FASTA file: '${PROTFILE}'"
	echo " - Output alignment file: '${OUTPUT}.aln'"

	mafft-linsi $PROTFILE > $OUTPUT".aln"  # *mod*


elif [ $ALNMETHOD == "mafft" ] ; then
	
	if ! command -v mafft &> /dev/null
	then	
		echo -e "Cannot find mafft binary, make sure to export the path";
		usage
		exit 1;
	fi

	echo " - Input protein FASTA file: '${PROTFILE}'"
	echo " - Output alignment file: '${OUTPUT}.aln'"

	mafft --auto $PROTFILE > $OUTPUT".aln"  # *mod*

else
	echo -e "\nERROR, no recognized MSA program was detected. Please specify 'mafft-linsi' or 'mafft' in -m option\n"
	exit 1;
fi


# Check if the aligment file is empty or not

if [ -s $OUTPUT".aln" ]
then
	echo -e "\nEnd"
	echo -e "S1.1 Finished successfully \n"
else
	echo "ERROR: Empty alignment file!!!"
	echo -e "S1.1 Finished with error \n"
	exit 1
fi

# Run phylogeny

echo "# | S1.2 Step |"
echo -e "# - Run phylogeny with ${TREEMETHOD}- #\n"
echo -e "Start\n"

if [ $TREEMETHOD == "FastTree" ] ; then

#	if [[ ! -f FastTree ]] ; then
	if ! command -v FastTree &> /dev/null
	then
		echo -e "Cannot find FastTree binary, make sure to export the path";
		usage
		exit 1;
	fi

	echo " - Input alignment file: '${OUTPUT}.aln'"
	echo " - Output treefile file: '${OUTPUT}.${TREEMETHOD}.nwk'"

	FastTree $OUTPUT".aln" > ${OUTPUT}.$TREEMETHOD".nwk"

elif [ $TREEMETHOD == "iqtree" ] ; then
	
#	if [[ ! -f iqtree2 ]] ; then
	if ! command -v iqtree2 &> /dev/null
	then		
		echo -e "Cannot find FastTree binary, make sure to export the path";
		usage
		exit 1;
	fi

	echo " - Input alignment file: '${OUTPUT}.aln'"
	echo " - Output treefile file: '${OUTPUT}.${TREEMETHOD}.nwk'"

	iqtree2 -s $OUTPUT".aln" -m MFP --prefix $OUTPUT
	mv $OUTPUT".treefile" ${OUTPUT}.$TREEMETHOD".nwk"

else
	echo -e "\nERROR, no recognized tree program was detected. Please specify 'FastTree' or 'iqtree' in -t option\n"
	exit 1;
fi

# Check if the newick file is empty or not

if [ -s ${OUTPUT}.$TREEMETHOD".nwk" ]
then
	echo -e "\nEnd"
	echo -e "S1.2 Finished successfully \n"
else
	echo "ERROR: Empty newick file!!!"
	echo -e "S1.2. Finished with error \n"
	exit 1
fi

# Get distances

#if [[ ! -f nw_distance ]] ; then
if ! command -v nw_distance &> /dev/null
then	
	echo -e "Cannot find nw_distance binary, make sure to export the path";
	usage
	exit 1;
fi

echo "# | S1.3 Step |"
echo -e "# - Evolutionary distances estimation - #\n"
echo -e "Start\n"

echo " - Input treefile file: '${OUTPUT}.$TREEMETHOD.nwk'"
echo " - Output distnace matrix: '${OUTPUT}.${TREEMETHOD}._distancematrix.tsv'"

nw_distance -n -m m $OUTPUT.$TREEMETHOD".nwk" > $OUTPUT.$TREEMETHOD"_distancematrix.tsv"

if [ -s $OUTPUT.$TREEMETHOD"_distancematrix.tsv" ]
then
	echo -e "\nEnd"
	echo -e "S1.3 Finished successfully \n"
else
	echo "ERROR: Empty *distancematrix.tsv file!!!"
	echo -e "S1.3. Finished with error \n"
	exit 1
fi

# echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

