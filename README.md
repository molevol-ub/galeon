# GALEON
#### A Comprehensive Bioinformatic Tool to Analyse and Visualise Gene Clusters in Complete Genomes

<div>
    <img src="https://github.com/molevol-ub/galeon/blob/main/GALEON_masterScripts/GaleonLogo.png" alt="Software Logo" width="200">
</div>

To facilitate the identification, analysis, and visualisation of physically clustered gene family genes within chromosome-level genomes, we introduce GALEON, a user-friendly bioinformatic tool. GALEON identifies gene clusters by studying the spatial distribution of pairwise physical distances among gene family members along with the genome-wide gene density. The pipeline also enables the simultaneous analysis and comparison of two gene families, and allows the exploration of the relationship between physical and evolutionary distances. This tool offers a novel approach for studying the origin and evolution of gene families.


GALEON documentation can be also be found in: http://www.ub.edu/softevol/galeon

### Version history

V1: Initial release

## 0. Contents

 1. Installation and prerequisites
 2. Input data
 3. Running GALEON
 4. Output
 5. Example data
 6. Citation
 7. Troubleshooting


## 1. Installation and prerequisites

GALEON is distributed as a set of several scripts that can be called from Galeon_masterScripts folder, but do not require any specific installation or compiilation step. However, the pipelines does require the following python modules, R packages as well as some additional software.

### 1.1. Python Packages
If conda is avaiilable, all these packages can be easily obtained by installing the Galeon conda environment, as well as some of the required software: mafft and newick_utils. Alternatively, you may install them separately using `pip`, consult the appropriate documentation for each of them.

```
argparse
ast
collections
copy
gc
itertools
matplotlib
numpy
operator
os
pandas
re
seaborn
shutil
string
scipy
subprocess
sys
time
```

### 1.2. R Packages

Make sure to have R installed, as well as two additional R packages: `rmarkdown` and `DT`, which are needed to create the final Report in HTML format. These packages can be installed as follows directly from the R terminal.


```
# 1-Open a terminal and run “R”
R

# 2-Install packages
>install.packages(“rmarkdown”)
>install.packages(“DT”)

# 3-Check that they can be loaded
>library(“rmarkdown”)
>library(“DT”)
```

### 1.3. Additional software

The following programs must be installed and available from command line: **mafft**, **bedtools**, **FastTree** and **iqtree2**. 

We provide a `bin` directory with binaries of `bedtools`, `FastTree` and `iqtree2`. If the Galeon conda environment will be created, `mafft` and `newick_utils` should be available upon environment activation. Alternatively, check the corresponding documentation for installation instructions.


**Notes**:
- **FastTree** - By default, GALEON uses FastTree to infer the phylogeny among gene family copies. http://www.microbesonline.org/fasttree/
- **Iqtree** - IqTree can be used instead of FastTree to reconstruct the gene family phylogeny: http://www.iqtree.org/

Tested software versions:

- Mafft v7.3.10
- bedtools v2.30.0
- FastTree v2.1.11
- iqtree2  v2.1.3
- python 3.11, 3.12
- R v4.1.2, v4.2.3


#### 1.4. Install GALEON

```
# 1-Download the software
git clone https://github.com/molevol-ub/galeon.git

# 2-Decompress the source code file (if needed)
cd galeon
tar -xf GALEON.v1.tar.gz 
cd GALEON_v1

# 3-Make the binaries executable
chmod +x GALEON_masterScripts/bin/*

# 4-Activate conda and install the Galeon conda environment
conda activate
conda create env -f GaleonEnv.yml

# 5-Activate the environment
conda activate Galeon

# 6-Run the configuration script
# this will add a header like this “#!/home/user/miniconda3/envs/Galeon/bin/python” to the python scripts
python Configure.py YOURPATH_to/GALEON_masterScripts

```

**Dependencies installation checkpoint**

Once all the packages have been installed, run the following command to check that all the dependencies are available and accessible. 

```
# 7-Enter to the GALEON_masterScripts directory and run the following script
cd GALEON_masterScripts
python Scripts/Check_installed_packages_and_PythonEnv.py 
```

If you encounter any errors related to the software (bedtools, mafft, iqtree2, FastTree), check the help message to add the path to your own installation

```
python Scripts/Check_installed_packages_and_PythonEnv.py -h
```

**Export GALEON to PATH**

GALEON scripts should be preferably added to PATH to have general access.

```
# 8-Export the path_to_GALEON to PATH
nano ~/.bashrc

# add this line: export PATH=YOURPATH_to/GALEON_v1:$PATH
# save and exit
# run
source ~/.bashrc

# 9-Check the accessibility to the Galeon control script
which GALEON_ControlScript.py
# now it should output: YOURPATH_to/GALEON_ControlScript.py
```

## 2. Input data
**Warning**: Please be careful while preparing the inputs, follow precisely the below instructions. Input file names and formats are mandatory. 

GALEON uses three types of files for each gene family:
- Annotation files: with the coordinates of the genes.
- Proteins (or MSA): in FASTA format. Needed for evolutionary distances computation.
- Chromosome/Scaffold sizes file. Needed for the summary plots and final report.

### 2.1. Annotation files format

- Input coordinate file format: *GFF3*, *BED1*, *BED2* (check the formats below)
- Input coordinates file name: **{FAMILYNAME}_fam.{FORMAT}**
  - Examples:
   - *GFF3* file name: GR_fam.gff3
   - *BED1* file name: GR_fam.bed1
   - *BED2* file name: GR_fam.bed2

All the input coordinate files MUST be provided in the same file format.

#### 2.1.1. GFF3 format
- 9 tab separated columns. Only five of them will be used: `scaffold`, `feature` `start`, `end` and `attribute`.
- Make sure that the "gene" is present in the `feature` column.
- Gene names are inferred from the `attribute` column.
- The header is shown to clarify the meaning of each field but it is not required.

| scaffold | source | feature | start | end | score | strand | frame | attribute |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| Scaffold_14804_HRSCAF_18385 | AnnotGFF | gene | 41841903 | 41843055 | 0.51 | - | . | ID=g10232;blastphmmer;annot;Pos:1-383 |
| Scaffold_14804_HRSCAF_18385 | AnnotGFF | mRNA | 41841903 | 41843055 | 0.51 | - | . | ID=g10232.t1;Parent=g10232;blastphmmer;annot;Pos:1-383 |
| Scaffold_14804_HRSCAF_18385 | AnnotGFF | CDS | 41841904 | 41843052 | 0.51 | - | 0 | ID=g10232.t1.CDS1;Parent=g10232.t1;blastphmmer;annot;Pos:1-383 |
| Scaffold_14804_HRSCAF_18385 | AnnotGFF | gene | 47268322 | 47268742 | 0.67 | + | . | ID=g10331;blastphmmer;annot;Pos:1-139 |
| Scaffold_14804_HRSCAF_18385 | AnnotGFF | mRNA | 47268322 | 47268742 | 0.67 | + | . | ID=g10331.t1;Parent=g10331;blastphmmer;annot;Pos:1-139 |
| Scaffold_14804_HRSCAF_18385 | AnnotGFF | CDS | 47268322 | 47268738 | 0.67 | + | 0 | ID=g10331.t1.CDS1;Parent=g10331.t1;blastphmmer;annot;Pos:1-139 |
| Scaffold_14804_HRSCAF_18385 | AnnotGFF | gene | 47277448 | 47277868 | 0.63 | + | . | ID=g10332;blastphmmer;annot;Pos:1-139 |
| Scaffold_14804_HRSCAF_18385 | AnnotGFF | mRNA | 47277448 | 47277868 | 0.63 | + | . | ID=g10332.t1;Parent=g10332;blastphmmer;annot;Pos:1-139 |
| Scaffold_14804_HRSCAF_18385 | AnnotGFF | CDS | 47277448 | 47277864 | 0.63 | + | 0 | ID=g10332.t1.CDS1;Parent=g10332.t1;blastphmmer;annot;Pos:1-139 |

#### 2.1.2. BED2 format
- 4 tab separated columns
- This provides directly the coordinates and the gene IDs.

| Scaffold ID | start | end | attribute |
| ------------- | ------------- | ------------- | ------------- |
| Scaffold_14804_HRSCAF_18385 | 41841903 | 41843055 | g10232 |
| Scaffold_14804_HRSCAF_18385 | 47268322 | 47268742 | g10331 |
| Scaffold_14804_HRSCAF_18385 | 47277448 | 47277868 | g10332 |


#### 2.1.3 - BED1 format
- 3 tab separated columns
- This format does not include gene IDs. Gene names will be given by the gene order. Consider this example annotation file name: `GR_fam.bed1`. The genes will be named as *GR_0*, *GR_1*, etc.


| Scaffold ID | start | end |
| ------------- | ------------- | ------------- |
| Scaffold_14804_HRSCAF_18385 | 41841903 | 41843055 |
| Scaffold_14804_HRSCAF_18385 | 47268322 | 47268742 |
| Scaffold_14804_HRSCAF_18385 | 47277448 | 47277868 |
| Scaffold_14804_HRSCAF_18385 | 51844347 | 51844998 |
| Scaffold_14804_HRSCAF_18385 | 52310537 | 52311098 |

**NOTE**: This format might be of use when there is some kind of problem related with the gene names format. Then, the user may run some tests to check whether the input gene family is organized in cluster.


### 2.2. Proteins and MSA files format
To compute the evolutionary distances, you will need to provide either the proteins of your gene family of interest or the corresponding MSA in FASTA format.

- Input proteins FASTA and MSA format: **{FAMILYNAME}_fam.{FORMAT}**
  - Examples:
   - *Proteins* file name: GR_fam.fasta
   - *MSA* file name: GR_fam.aln

**NOTES**:

- Protein names MUST coincide with the gene name of the input GFF3 or BED2 file.
- If raw protein data is provided, *mafft* will align them.
- If pre-computed MSA data is provided, it will be used directly to run *FastTree* or *iqtree* and get the evo. distances.
- BED1 format does not contain the gene name information, so BED2 format should be used instead.
 

### 2.3. Chromosome/Scaffold sizes file
This is used mainly as a guide to filter the output results and summarise the findings focusing on the main scaffolds (those corresponding to chromosomes) or a subset of scaffolds of choice (for example: the ten largest scaffolds or a list of scaffolds of interest).

- Input file name: **ChrSizes.txt**
- 3 tab separated columns

| Scaffold ID | Length (in bp) | Scaffold associated name |
| ------------- | ------------- | ------------- |
| Scaffold_15362_HRSCAF_19823 | 317950935 | ChrX |
| Scaffold_14804_HRSCAF_18385 | 177171321 | Chr1 |
| Scaffold_14178_HRSCAF_16784 | 176727214 | Chr2 |



## 3. Running GALEON

### 3.1. Estimate g parameter (*mode: gestimate*)
In this mode, the pipeline estimates the expected number of genes found in a number of bases, as well as the number of genes expected across the g input values and the probability of finding 2 or more genes in a window of g size (i.e.: 100 Kb), which would be considered as a cluster in the following analyses (**Section 3.2.**).

Run the following command to estimate the g parameter based on the inputs. No input files are required here.

- `-n NUM` Gene family size (number of genes)
- `-s NUM` Genome size (in Mb units)
- `-g NUM` g value (in Kb units): more that one can be tested
- `-outdir DIRNAME` Output directory, set by default to *g_estimation_Results_Directory*.

**Commands**
```
# Run using one g value
GALEON_ControlScript.py gestimate -n 134 -s 1354 -g 100

# Test several g values
GALEON_ControlScript.py gestimate -n 134 -s 1354 -g 150,200,300,400
```

**Help message**
```
# Run using one g value
GALEON_ControlScript.py gestimate -h
```

**Results**

Estimation table: `g_estimation.table.txt`

Output table that contains:

column: "g value" - input g value

column: "Exp. 1 gene each X Mb" - 1 gene, expected to be found each "X" Mb

column: "Exp. genes / Mb" - # of genes, expected to be found each Mb

column: "Exp. Genes / Kb" - # of genes, expected to be found each 1 Kb

column: "Exp. Genes/g value" - # of genes, expected to be found each "g" Kb

column: "P(X>=2) / g value" - The probability of finding by chance two (or more) genes in a "g" kb stretch

column: "Poisson's λ"

Log files: 

Logs_gestimate_mode/gestimation.out 

Logs_gestimate_mode/gestimation.err



### 3.2. Gene cluster identification (*mode: clusterfinder*)




```
python ControlScript.py clusterfinder -a GFFs/ -g 100
```

*Continue manual*


## 3. Output 

Galeon will provide a summary output and plots in markdown... *fill*



## 4. Example dataset


An example to run GALEON can be found in Example_data folder. 

To run GALEON with the example dataset, the following commands should be used:

```
commands?
```

## 9. Citation


Vadim Pisarenco, Joel Vizueta, Julio Rozas. GALEON . Submitted. 2023.



## 10. Troubleshooting


If you find any error, please create an issue in Github specifying the error and all details as possible.


