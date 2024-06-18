# GALEON
#### A Comprehensive Bioinformatic Tool to Analyse and Visualise Gene Clusters in Complete Genomes

<div>
    <img src="https://github.com/molevol-ub/galeon/blob/main/GALEON_masterScripts/GaleonLogo.png" alt="Software Logo" width="200">
</div>


To facilitate the identification, analysis, and visualisation of physically clustered gene family genes within chromosome-level genomes, we introduce GALEON, a user-friendly bioinformatic tool. GALEON identifies gene clusters by studying the spatial distribution of pairwise physical distances among gene family members along with the genome-wide gene density. The pipeline also enables the simultaneous analysis and comparison of two gene families, and allows the exploration of the relationship between physical and evolutionary distances. This tool offers a novel approach for studying the origin and evolution of gene families.


GALEON documentation can be also be found in: http://www.ub.edu/softevol/galeon

### Version history

V1: Initial release

## Contents

 1. Installation and prerequisites
 2. Input data
 3. Running GALEON
 4. Example data
 5. Citation
 6. Troubleshooting

## 1. Installation and prerequisites

GALEON is distributed as a set of scripts that can be called from Galeon_masterScripts folder, but do not require any specific installation or compilation step. However, the pipeline does require several python modules, R packages, as well as external software. All of them are listed in Section 1.2, 1.3.

It is highly recommened to install the conda Galeon environment that provides all of the required python packages as well as some of the external programs, specifically pandoc and mafft.

### 1.1. Install GALEON

```
# 1-Download the software
git clone https://github.com/molevol-ub/galeon.git
cd galeon

# 2-Make the binaries executable
chmod +x GALEON_masterScripts/bin/*

# 3-Activate conda and install the Galeon conda environment
conda activate
conda env create -f GaleonEnv.yml

# 4-Activate the environment
conda activate Galeon

# 5-Run the configuration script
# this will add a header like this “#!/home/user/miniconda3/envs/Galeon/bin/python” to the python scripts
python Configure.py YOURPATH_to/GALEON_masterScripts

```

**Dependencies installation checkpoint**

Once all the packages have been installed, run the following command to check that all the dependencies are available and accessible. 

```
# 6-Enter to the GALEON_masterScripts directory and run the following script
cd GALEON_masterScripts
python Scripts/Check_installed_packages_and_PythonEnv.py 
```

If you encounter any errors related to the software (bedtools, mafft, iqtree2, FastTree), check the help message to add the path to your own installation. 

**In addition, note** that `R` and two R packages, `rmarkdown` and `DT`, need to be installed (see 1.3 to install them). 

```
python Scripts/Check_installed_packages_and_PythonEnv.py -h
```

**Export GALEON to PATH**

GALEON scripts should be preferably added to PATH to have general access.

```
# 7-Export the path_to_GALEON to PATH
nano ~/.bashrc

# add this line: export PATH=YOURPATH_to/GALEON_masterScripts:$PATH
# save and exit
# run
source ~/.bashrc

# 8-Check the accessibility to the Galeon control script
which GALEON_ControlScript.py
# now it should output: YOURPATH_to/GALEON_ControlScript.py
```

### 1.2. Python Packages
It is highly recommended to use conda since all the python packages will be easily installed with the Galeon conda environment (see Section 1.1), as well as some of the required software: mafft and newick_utils. Alternatively, you may install them separately using `pip`, consult the appropriate documentation for each of them.

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

### 1.3. R Packages

Make sure to have R installed, as well as two additional R packages: `rmarkdown` and `DT`, which are needed to create the final Report in HTML format. These packages can be installed as follows directly from the R terminal.


```
# 1-Open a terminal and run “R”
R

# 2-Install packages
>install.packages("rmarkdown")
>install.packages("DT")

# 3-Check that they can be loaded
>library("rmarkdown")
>library("DT")
```

### 1.4. Additional software

The following programs must be installed and available from command line: **pandoc**, **mafft**, **bedtools**, **FastTree** and **iqtree2**. 

We provide a `bin` directory with binaries of `bedtools`, `FastTree` and `iqtree2`. If the Galeon conda environment is created, `pandoc`, `mafft` and `newick_utils` should be available upon environment activation. 

**Notes**:
- **Iqtree** - To reconstruct the gene family phylogeny, [IqTree](http://www.iqtree.org/) can be run using with default presets or using '--fast' option in order to perform a fast tree search (resembling the FastTree method).
    - By default, GALEON runs iqtree with '--fast' option (`-t iqtree-fast`)
    - For more accurate tree search run GALEON changing the `-t` parameter (`-t iqtree`)
- **FastTree** - By default, GALEON uses FastTree to infer the phylogeny among gene family copies. http://www.microbesonline.org/fasttree/

Alternatively, check the corresponding documentation for installation instructions.

Tested software versions:

- Mafft v7.3.10
- bedtools v2.30.0
- FastTree v2.1.11
- iqtree2 v2.1.3
- python 3.11, 3.12
- perl v5.32.1
- R v4.1.2, v4.2.3
- pandoc v2.9.2

**Known issues**
- pandoc v3.1.13 works but the resulting report shows some non-optimal visualization. Conda Galeon environment uses pandoc v2.9.2.

## 2. Input data
**Warning**: Please be careful while preparing the inputs, we recommend to carefully read and follow the below instructions. Input file name structure and formats are mandatory.

GALEON uses three types of files for each gene family:
- 2.1 - Annotation files containing the coordinates of the genes.
- 2.2 - Proteins (or MSA) in FASTA format. Needed for evolutionary distances computation.
- 2.3 - Chromosome/Scaffold sizes file. Needed for the summary plots and final report.

### 2.1. Annotation files

- Input coordinate file format: *GFF3*, *BED1*, *BED2* (check the formats below)
- Input coordinates file name: **{FAMILYNAME}_fam.{FORMAT}**
  - Examples:
   - *GFF3* file name: GR_fam.gff3
   - *BED1* file name: GR_fam.bed1
   - *BED2* file name: GR_fam.bed2

All the input coordinate files MUST be provided in the same file format.

#### 2.1.1. GFF3 format
- It contains 9 tab separated columns, from which the following will be used in Galeon: `scaffold`, `feature` `start`, `end` and `attribute`.
- Make sure that "gene" is present in the `feature` column.
- Gene names are inferred from the `attribute` column (ID).
- The header is shown to clarify the meaning of each field but it is not required.

| scaffold | source | feature | start | end | score | strand | frame | attribute |
| :-------------: | :-------------: | :-------------: | :-------------: | :-------------: | :-------------: | :-------------: | :-------------: | ------------- |
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
- A file containing 4 tab separated columns.
- It provides directly the coordinates and the gene IDs.

| Scaffold ID | start | end | attribute |
| :-------------: | :-------------: | :-------------: | :-------------: |
| Scaffold_14804_HRSCAF_18385 | 41841903 | 41843055 | g10232 |
| Scaffold_14804_HRSCAF_18385 | 47268322 | 47268742 | g10331 |
| Scaffold_14804_HRSCAF_18385 | 47277448 | 47277868 | g10332 |


#### 2.1.3. BED1 format
- File with 3 tab separated columns.
- This format does not include gene IDs. Gene names will be given by the gene order. Consider this example annotation file name: `GR_fam.bed1`. The genes will be named as *GR_0*, *GR_1*, etc.


| Scaffold ID | start | end |
| :-------------: | :-------------: | :-------------: |
| Scaffold_14804_HRSCAF_18385 | 41841903 | 41843055 |
| Scaffold_14804_HRSCAF_18385 | 47268322 | 47268742 |
| Scaffold_14804_HRSCAF_18385 | 47277448 | 47277868 |
| Scaffold_14804_HRSCAF_18385 | 51844347 | 51844998 |
| Scaffold_14804_HRSCAF_18385 | 52310537 | 52311098 |

**NOTE**: This format might be useful when there is some kind of problem related with the gene names format. Then, the user may run some tests to check whether the input gene family is organized in cluster.


### 2.2. Proteins or MSA file
To compute the evolutionary distances, you will need to provide either the proteins of your gene family of interest or the corresponding MSA in FASTA format.

- Input proteins FASTA and MSA format: **{FAMILYNAME}_fam.{FORMAT}**
  - Examples:
   - *Proteins* file name: GR_fam.fasta
   - *MSA* file name: GR_fam.aln

**NOTES**:

- Protein names MUST coincide with the gene name of the input GFF3 or BED2 file.
- If raw protein data is provided, our pipeline will use *mafft* to align them.
- If pre-computed MSA data is provided, it will be used directly to run *FastTree* or *iqtree* and get the evolutionary distances.
- BED1 format does not contain the gene name information, so BED2 format should be used instead.
 

### 2.3. Chromosome/Scaffold size file
This is used mainly as a guide to filter the output results and summarise the findings focusing on the main scaffolds (those corresponding to chromosomes) or a subset of scaffolds of choice (for example: the ten largest scaffolds or a list of scaffolds of interest).

- Input file name: **ChrSizes.txt**
- 3 tab separated columns

| Scaffold ID | Length (in bp) | Scaffold associated name |
| :-------------: | :-------------: | :-------------: |
| Scaffold_15362_HRSCAF_19823 | 317950935 | ChrX |
| Scaffold_14804_HRSCAF_18385 | 177171321 | Chr1 |
| Scaffold_14178_HRSCAF_16784 | 176727214 | Chr2 |

If you don't have this file, you can generate by running the following script: `Get_scaffold_length.pl` located in `GALEON_masterScripts/Scripts`. It takes as an input your genome file in FASTA format. 

- Input file: YOURgenome.fasta
- Output file: YOURgenome_scaffold_length_sorted.txt

Command example:

```
perl GALEON_masterScripts/Scripts/Get_scaffold_length.pl YOURgenome.fasta
```

NOTE: The output table will have the above-described 3-column format. Optionally, you can rename the third column by replacing the *Scaffold IDs* with *Scaffold associated names*. Check the [example file](https://github.com/molevol-ub/galeon/blob/main/Example_data/Test_2/ChrSizes.txt).

<br>

## 3. Running GALEON

### 3.1. Estimate *g* parameter (*mode: gestimate*)
In this mode, the pipeline estimates the expected number of genes found in a number of bases, as well as the number of genes expected across the *g* input values and the probability of finding 2 or more genes in a window of g size (i.e.: 100 Kb), which would be considered as a cluster in the following analyses (**Section 3.2.**).

Run the following command to estimate the g parameter based on the inputs. No input files are required here.

- `-n NUM` Gene family size (number of genes)
- `-s NUM` Genome size (in Mb units)
- `-g NUM` g value (in Kb units): more than one can be tested
- `-outdir DIRNAME` Output directory, set by default to *g_estimation_Results_Directory*.

**Help message**

```
GALEON_ControlScript.py gestimate -h
```

**Commands**

```
# Run using one g value
GALEON_ControlScript.py gestimate -n 134 -s 1354 -g 100

# Test several g values
GALEON_ControlScript.py gestimate -n 134 -s 1354 -g 150,200,300,400
```

**Output**

```
# Output table
g_estimation_Results_Directory/
└── g_estimation.table.txt

# Logs and error messages
Logs_gestimate_mode/
├── gestimation.err
└── gestimation.out
```
- Output table: `g_estimation.table.txt`
The table that contains:

1. column: **g value** - input g value
2. column: **Mb per gene family member** - 1 gene, expected to be found each "X" Mb
3. column: **Expected genes per Mb** - # of genes, expected to be found each Mb
4. column: **Exp. genes per g value** - # of genes, expected to be found each 1 Kb
5. column: **Exp. Genes/g value** - # of genes, expected to be found each "g" Kb
6. column: **P(X>=2) in g kb** - The probability of finding by chance two (or more) genes in a "g" kb stretch
7. column: **Poisson's lambda**

<br>

### 3.2. Gene cluster identification (*mode: clusterfinder*)
In this mode the pipeline analyzes one (or several) gene families' data to identify clusters of genes in the genome. 

**Help message**
```
GALEON_ControlScript.py gestimate -h
```

#### 3.2.1. Single family analysis using physical distances

In this case, coordinates files are going to be analyzed to get pairwise distances between genes and arrange them in a distance matrix. This matrix will then be scanned to identify gene clusters. Finally, the distance matrix will be displayed as a heatmap with all the identified clusters (if any) represented by black square shapes.

**Inputs**

- More than one families can be analyzed at the same time

```
# How your input annotation directory "GFFs/" should look
└── GFFs
    ├── GR_fam.gff3
    └── IR_fam.gff3
```

**Commands**

Follow the instruction to run the analysis, generate the plots and create a final portable HTML report which will provide an overview of all the obtained results at a glance.

**Step 1)** Find clusters, independently for each input gene family using the coordinates files.

```
# Simplest command to run Galeon
GALEON_ControlScript.py clusterfinder -a GFFs/ -e disabled
```

- `-a DIRNAME` Input annotation directory with the coordinate files `-a GFFs`)
- `-g NUM` The g value is set to 100 Kb by default (`-g 100`)
- `-o DIRNAME` Output directory name, set by default to `-o clusterfinder_Results_Directory`

NOTE: The use of evolutionary distances is disabled here (`-e disabled`).

**Step 2)** Generate summary plots and tables for each input gene family

```
# Generate summary files for the GR family
GALEON_SummaryFiles.py -fam GR -clust clusterfinder_Results_Directory/ -coords GFFs -ssize ChrSizes.txt -sfilter 7

# Generate summary files for the IR family
GALEON_SummaryFiles.py -fam IR -clust clusterfinder_Results_Directory/ -coords GFFs -ssize ChrSizes.txt -sfilter 7
```

- `-fam FAMILYNAME` Your gene family name. Remember how is the format (**{FAMILYNAME}_fam.{FORMAT}**). For ex: "GR_fam.fasta" => "GR"
- `-clust DIRNAME` Input Galeon results directory. For ex: `clusterfinder_Results_Directory`
- `-coords DIRNAME` Annotation directory with the coordinate files, same as in `-a GFFs`
- `-ssize FILE` Chromosome/Scaffold size file.
- `-sfilter ALL|NUM|FILE` The summary plots will represent the results for a "NUM" number of largest scaffolds; a list of scaffolds of interest provided as a single column in an input "FILE"; "ALL" the scaffolds (often too many, the resulting summary plots might not informative).

**Step 3)** Generate a final HTML report
```
# Generate the final HTML report
GALEON_Report.py -clust clusterfinder_Results_Directory/ -ssize ChrSizes.txt -echo False
```

- `-clust DIRNAME` Input Galeon results directory. For ex: `clusterfinder_Results_Directory`
- `-ssize FILE` Chromosome/Scaffold size file.
- `-echo True|False` If `True` complete the paths to the files and plots are shown in the report. 

<br>

**Output**

GALEON will generate a portable HTML report, one for each family and per tested g value. It will contain all the generated tables and reports, making it easy to quickly access all the results. The report includes a "HELP" tab with useful information for interpreting the results.

Galeon results directory content

- **PhysicalDist_Matrices** directory
    - This directory containts pairwise physical distance matrices (`*matrix`) between genes of a given Gene Family, in each scaffold.
    - Also, physical distance heatmaps are present in svg and pdf format.
- **Plots** directory
    - This directory contains summary plots with cluster size distribution at genome level (dir.: `SummaryPlots_100.0g`) and individual scaffold level (dir.: `IndividualPlots_100.0g`)
    - Also, it contains some summary tables.
        - `X_family_ClusterSizes.table.100.0g.tsv` - This table contains a more explicit information about the size of each cluster in each scaffold.
        - `X_family_GeneLocation.table.100.0g.tsv` - It contains information detailed at gene level, specifically, gene coordinates and membership (singleton or clustered).
        - `X_family_GeneOrganizationGenomeSummary.table.100.0g.tsv` - Table with the number of genes organized in clusters (Clustered) or not (Singleton category) at the genome level.
        - `X_family_GeneOrganizationSummary.table.100.0g.tsv` - This table informs about how many genes are organized in clusters (Clustered) or are physically isolated (Singleton category) for each scaffold.

- **Reports** directory
    - This directory contains html reports, one for each family and g value.


Galeon results directory, example tree-like representation for two gene families "GR_fam" and "IR_fam":

```
clusterfinder_Results_Directory/
├── PhysicalDist_Matrices
│   ├── GR_fam.gff3.temp_matrices # Physical distance matrices, *matrix (in bp units)
│   ├── IR_fam.gff3.temp_matrices 
│   ├── GR_fam.gff3.temp_matrices_100.0g # Heatmaps in svg and pdf format
│   └── IR_fam.gff3.temp_matrices_100.0g
|
├── Plots # Contains summary tables and plots for each input family
│   ├── GR_fam
│   │   ├── GR_family_ClusterSizes.table.100.0g.tsv
│   │   ├── GR_family_GeneLocation.table.100.0g.tsv
│   │   ├── GR_family_GeneOrganizationGenomeSummary.table.100.0g.tsv
│   │   ├── GR_family_GeneOrganizationSummary.table.100.0g.tsv
│   │   |
│   │   ├── IndividualPlots_100.0g 
│   │   └── SummaryPlots_100.0g
│   └── IR_fam
│   │   ├── IR_family_ClusterSizes.table.100.0g.tsv
│   │   ├── IR_family_GeneLocation.table.100.0g.tsv
│   │   ├── IR_family_GeneOrganizationGenomeSummary.table.100.0g.tsv
│   │   ├── IR_family_GeneOrganizationSummary.table.100.0g.tsv
│   │   |
│   │   ├── IndividualPlots_100.0g
│   └── └── SummaryPlots_100.0g
|
└── Reports # One report for each family and tested g value.
    ├── GR_fam_100.0g_Report.html
    └── IR_fam_100.0g_Report.html
```

<br>

#### 3.2.2. Single family analysis using physical and evolutionary distances

Coordinates files will be processed as described in Section 3.2.1 to obtain the matrices and identify the clusters. However, in addition to coordinate files, proteins will be included to compute evolutionary distances. These distances will then be merged with the physical distance matrix by replacing the upper semi-matrix values. This "merged" matrix will be displayed as a heatmap, with all the identified clusters (if any) represented by black square shapes.

**Inputs**

- More than one families can be analyzed at the same time
- Provide the protein data in FASTA format, either as raw proteins or multiple sequence alignment (MSA). Plesae, use only ONE format.

```
# How your input annotation directory "GFFs/" and protein directory "Proteins/" should look
├── GFFs
│   ├── GR_fam.gff3
│   └── IR_fam.gff3
└── Proteins
    ├── GR_fam.fasta # or GR_fam.aln if MSA is provided
    └── IR_fam.fasta # or IR_fam.aln
```

**Commands**

Follow the instruction to run the analysis, generate the plots and create a final portable HTML report which will provide an overview of all the obtained results at a glance.

**Step 1)** Find clusters using the input coordinates and protein files 

```
# Simplest command to run Galeon
# Run this...
GALEON_ControlScript.py clusterfinder -a GFFs/ -e enabled -p Proteins/

# ...or this if the MSA files are already present in the "Proteins/" directory for each of the gene families of interest
GALEON_ControlScript.py clusterfinder -a GFFs/ -e enabled -p Proteins -pm True
```

- `-a DIRNAME` Input annotation directory with the coordinate files `-a GFFs`)
- `-g NUM` The g value is set to 100 Kb by default (`-g 100`)
- `-e enabled|disabled` Enabled the use of proteins (`-e enabled`)
- `-p DIRNAME` Input directory with Proteins or MSA files (`-p Proteins`)
- `-o DIRNAME` Output directory name, set by default to `-o clusterfinder_Results_Directory`
- `-feat gene|mRNA` ID used in the protein sequences that matches the gff3. The gene id is read by default (`-feat mRNA`) 

**NOTE:** Remember that the Protein and Gene IDs must be equal. For example, let's consider a GFF3 file with a "gene" named "ABC" and "mRNA" named "ABC.t1". 

| scaffold | source | feature | start | end | score | strand | frame | attribute |
| :-------------: | :-------------: | :-------------: | :-------------: | :-------------: | :-------------: | :-------------: | :-------------: | ------------- |
| Scaffold1 | AnnotGFF | gene | 100 | 1000 | . | - | . | ID=ABC;annot;Pos:1-409; |
| Scaffold1 | AnnotGFF | mRNA | 100 | 1000 | . | - | . | ID=ABC.t1;Parent=ABC;annot;Pos:1-409; |

- If your protein ID matches the "mRNA" ID, use the above commands. 
- But, if the protein ID matches the "gene" ID, then you will need to set the paramter `-feat gene` (by default, it is set to `-feat mRNA`).

```
# (modified) Simplest command to run Galeon
# Run this...
GALEON_ControlScript.py clusterfinder -a GFFs/ -e enabled -p Proteins/ -feat gene

# (modified) ...or this if the MSA files are already present in the "Proteins/" directory for each of the gene families of interest
GALEON_ControlScript.py clusterfinder -a GFFs/ -e enabled -p Proteins -pm True -feat gene
```

**Step 2)** Get evolutionary statistics (Cst) and perform the Mann-Whitney test
```
# Get evolutionary statistics (Cst) and perform the Mann-Whitney test
GALEON_GetEvoStats.py -clust clusterfinder_Results_Directory/ -prot Proteins/ -coords GFFs
```

**Step 3)** Generate summary plots, tables and the HTML report

```
# Generate summary files for the GR family
GALEON_SummaryFiles.py -fam GR -clust clusterfinder_Results_Directory/ -coords GFFs -ssize ChrSizes.txt -sfilter 7

# Generate summary files for the IR family
GALEON_SummaryFiles.py -fam IR -clust clusterfinder_Results_Directory/ -coords GFFs -ssize ChrSizes.txt -sfilter 7

# Generate the final HTML report
GALEON_Report.py -clust clusterfinder_Results_Directory/ -ssize ChrSizes.txt -echo False
```

<br>


**Output**

Galeon results directory content
- **MergedDistances_Dataframes** directory.
    - This directory containts pairwise physical and evolutionary distance files (`*matrix` and `*distance.tsv`) between genes of a given Gene Family, in each scaffold.
    - Also, it contains physical and evolutionary distance heatmaps and scatterplots in svg and pdf format at a genome level (representing the entire gene family) as well as for each scaffold (showing each one separately).
- **PhysicalDist_Matrices** directory. Same as in Section 3.2.1.
- **Plots** directory. Same as in Section 3.2.1.
- **Reports** directory. Same as in Section 3.2.1.

Galeon results directory, example tree-like representation for two gene families "GR_fam" and "IR_fam":

```
clusterfinder_Results_Directory/
│
├── MannWhitney_StatisticsResults
│    ├── GR_fam_GlobalStats_value.100.0g.txt
│    └── GR_fam_MannWhitney.results.brief.100.0g.tsv
│
├── MergedDistances_Dataframes
│    ├── GR_fam.IntermediateFiles/
│    ├── GR_fam.merged.matrices/ # Physical + Evolutionary distance matrices
│    ├── GR_fam.plots_100.0g/ # Physical + Evolutionary distance Heatmaps and Scatterplots in svg and pdf format
│    │
│    ├── IR_fam.IntermediateFiles/
│    ├── IR_fam.merged.matrices/
│    ├── IR_fam.plots_100.0g/
│    │
│    ├── GR_fam.GlobScatterPlot_100.0g.pdf # Physical vs Evo. distance scatter plot considering at genome level, that is, considering all the genes of the input gene family
│    ├── GR_fam.GlobScatterPlot_100.0g.svg
│    ├── IR_fam.GlobScatterPlot_100.0g.pdf
│    └── IR_fam.GlobScatterPlot_100.0g.svg
│
├── PhysicalDist_Matrices
│    ├── GR_fam.gff3.temp_matrices/  # Physical distance matrices, *matrix (in bp units)
│    ├── GR_fam.gff3.temp_matrices_100.0g/ # Heatmaps in svg and pdf format
│    │
│    ├── IR_fam.gff3.temp_matrices/
│    └── IR_fam.gff3.temp_matrices_100.0g/
│
├── Plots # same content as in "Section 3.2.1"
│    ├── GR_fam
│    │   ├── IndividualPlots_100.0g/
│    │   └── SummaryPlots_100.0g/
│    └── IR_fam
│        ├── IndividualPlots_100.0g/
│        └── SummaryPlots_100.0g/
│
└── Reports/ # same content as in "Section 3.2.1"
    ├── GR_fam_100.0g_Report.html
    └── IR_fam_100.0g_Report.html

```

#### 3.2.3. Joint family analysis
Find clusters between two input families using the coordinates from the input files. Note that only two families can be analyzed at once, and protein sequences cannot be used in this mode.

**Inputs**

- Add the two gene family gff3 files that will be analyzed.

```
# How your input annotation directory "GFFs/" should look
└── GFFs
    ├── GR_fam.gff3
    └── IR_fam.gff3
```

**Commands**

**Step 1)** Find clusters using the input coordinates

- Remember that in this analysis, only two families must be present in the input annotation directory (`-a GFFs`).

```
# Simplest command to run Galeon
GALEON_ControlScript.py clusterfinder -a GFFs/ -e disabled -F BetweenFamilies
```

- `-a DIRNAME` Input annotation directory with the coordinate files `-a GFFs`)
- `-e enabled|disabled` Disable the use of proteins (`-e disabled`)
- `-F WithinFamilies|BetweenFamilies` Perform a separate analysis for each family or a joint analysis for two input families (`-F BetweenFamilies`)

**Step 2-3)** Generate summary plots, tables and the HTML report

```
# Generate summary files for the GR, IR families and merged data of both
GALEON_SummaryFiles.py -fam merged -clust clusterfinder_Results_Directory/ -coords GFFs/merged_dir/ -ssize ChrSizes.txt -sfilter 7

# Generate the final HTML report
GALEON_Report.py -clust clusterfinder_Results_Directory/ -ssize ChrSizes.txt -echo False
```

**Note** how the `-coord` parameter is specified here (`-coord GFFs/merged_dir`), it is a bit different because you must add "merged_dir" at the end of the command.

Galeon results directory content

- **PhysicalDist_Matrices** directory
    - This directory containts pairwise physical distance matrices (`*matrix`) between genes of a given Gene Family, in each scaffold.
    - Also, physical distance heatmaps are present in svg and pdf format.
- **Plots** directory
    - This directory contains summary plots with cluster size distribution at genome level (dir.: `SummaryPlots_100.0g`) and individual scaffold level (dir.: `IndividualPlots_100.0g`)
    - Since now we are interested in the joint analysis of two families (let's name them X and Y), the output summary plots are presented in 4 subdirectories:
        - `X` - plots for clusters of X gene family
        - `Y` - plots for clusters of Y gene family
        - `X.Y` - plots for clusters of X and Y gene family members (when members both families belong to the same cluster).
        - `merged` - plots for clusters of X and Y gene family considered as one.
    - Also, it contains some summary tables.
        - `X_family_ClusterSizes.table.100.0g.tsv` - This tables contains a more explicit information about the size of each cluster in each scaffold.
        - `merged_family_GeneLocation.table.100.0g.tsv` - This tables contains information detailed at gene level of both families, specifically, gene coordinates and membership (singleton or clustered).
        - `X_family_GeneOrganizationGenomeSummary.table.100.0g.tsv` - This tables informs about how many genes are organized in clusters (Clustered category) or are not (Singleton category) at genome level.
        - `X_family_GeneOrganizationSummary.table.100.0g.tsv` - This tables informs about how many genes are organized in clusters (Clustered category) or are not (Singleton category) in each scaffold.
        - `X.Y_family_GeneOrganizationSummary.table.100.0g.tsv` - This tables informs about how many genes of both families are jointly organized in clusters (Clustered category) or are not (Singleton category) in each scaffold. This file appears only when there are clusters formed by members of both families. Otherwise, it is absent.

- **Reports** directory
    - This directory contains html reports, one for each family and g value.

Galeon results directory, example tree-like representation for two gene families "GR_fam" and "IR_fam":

```
clusterfinder_Results_Directory/
├── PhysicalDist_Matrices
│   ├── merged_fam.gff3.temp_matrices # contains matrices *matrix
│   └── merged_fam.gff3.temp_matrices_100.0g # physical distance plots in svg and pdf format
│
├── Plots
│   └── merged_fam # summary tables and plots
│       ├── GR_family_ClusterSizes.table.100.0g.tsv
│       ├── GR_family_GeneOrganizationGenomeSummary.table.100.0g.tsv
│       ├── GR_family_GeneOrganizationSummary.table.100.0g.tsv
│       │
│       ├── IR_family_ClusterSizes.table.100.0g.tsv
│       ├── IR_family_GeneOrganizationGenomeSummary.table.100.0g.tsv
│       ├── IR_family_GeneOrganizationSummary.table.100.0g.tsv
│       │
│       ├── GR.IR_family_GeneOrganizationSummary.table.100.0g.tsv
│       │
│       ├── merged_family_GeneLocation.table.100.0g.tsv
│       │
│       ├── IndividualPlots_100.0g
│       │   ├── GR # Gr family clusters' size distribution
│       │   ├── GR.IR # Two family clusters size distribution
│       │   ├── IR # Gr family clusters' size distribution
│       │   └── merged # All clusters' size distribution
│       └── SummaryPlots_100.0g
│           ├── GR
│           ├── GR.IR
│           ├── IR
│           └── merged
│           
└── Reports
    ├── merged_fam_100.0g_Report.html
    └── merged_fam_100.0g_Report.Rmd
```

<br>


## 4. Example dataset

Several test datasets are available in [Example_data](https://github.com/molevol-ub/galeon/tree/main/Example_data) directory. You can enter to each *test* directory and run the commands described below. The output files for each example are also provided in the file [Solved_examples.tar.gz](https://github.com/molevol-ub/galeon/blob/main/GALEON_compressed_file/Solved_examples.tar.gz)
- Remember to have the Galeon conda environment activated (or alternativelly all the dependencies installed) before running the commands.


### 4.1. Test 1. Estimate g parameter

*No input files are required here.
- Gene family size (number of genes): `-n 134`
- Genome size (in Mb units): `-s 1354`
- g value (in Kb units): one `-g 100` or several values `-g 100,200,300`

```
cd Test_1

# Run using one g value
GALEON_ControlScript.py gestimate -n 134 -s 1354 -g 100

# Test several g values
GALEON_ControlScript.py gestimate -n 134 -s 1354 -g 150,200,300,400
```

- Check the output directory "g_estimation_Results_Directory"

### 4.2. Test 2. Single family analysis using physical distances

- Input directory with annotation files: `-a GFFs`
- Input g value: `-g 100,200`
- Output directory: `-outdir 2_OneFam_PhysDistOnly_GFF3`
- Chromosome/Scaffold size file: `-ssize ChrSizes.txt`
- Summarize the results for the first 7 largest scaffolds: `-sfilter 7`
- Write the path to files in the final report?: `-echo False`

```
cd Test_2

# Run Galeon
GALEON_ControlScript.py clusterfinder -a GFFs/ -g 100,200 -e disabled -outdir 2_OneFam_PhysDistOnly_GFF3

# Generate summary files and tables for GR family
GALEON_SummaryFiles.py -fam GR -clust 2_OneFam_PhysDistOnly_GFF3/ -coords GFFs -ssize ChrSizes.txt -sfilter 7

# Generate summary files and tables for IR family
GALEON_SummaryFiles.py -fam IR -clust 2_OneFam_PhysDistOnly_GFF3/ -coords GFFs -ssize ChrSizes.txt -sfilter 7

# Create a summary report
GALEON_Report.py -clust 2_OneFam_PhysDistOnly_GFF3/ -ssize ChrSizes.txt -echo False
```

### 4.3. Test 3. Single family analysis using physical and evolutionary distances (using unaligned protein sequences)
- Input directory with annotation files: `-a GFFs`
- Input g value: `-g 100,200`
- Enable the usage of Proteins: `-e enabled`
- Protein directory with fasta files: `-p Proteins`
- Output directory: `-outdir 3_OneFam_PhysEvoDistances_GFF3`
- Square frame color: `-f orange`
- Chromosome/Scaffold size file: `-ssize ChrSizes.txt`
- Summarize the results for the first 7 largest scaffolds: `-sfilter 7`
- Write the path to files in the final report?: `-echo False`

There are three directories, one for each option to compute the evolutionary distance:
- Test_3_iqtree-fast
- Test_3_iqtree
- Test_3_FastTree

```
cd Test_3_iqtree-fast/

# Run Galeon
GALEON_ControlScript.py clusterfinder -a GFFs/ -g 100,300 -e enabled -p Proteins/ -outdir 3_OneFam_PhysEvoDistances_GFF3/ -f orange

# Get evolutionary statistics (Cst) and perform the Mann-Whitney test
GALEON_GetEvoStats.py -clust 3_OneFam_PhysEvoDistances_GFF3/ -prot Proteins/ -coords GFFs

# Generate summary files and tables for GR family
GALEON_SummaryFiles.py -fam GR -clust 3_OneFam_PhysEvoDistances_GFF3/ -coords GFFs -ssize ChrSizes.txt -sfilter 7

# Generate summary files and tables for IR family
GALEON_SummaryFiles.py -fam IR -clust 3_OneFam_PhysEvoDistances_GFF3/ -coords GFFs -ssize ChrSizes.txt -sfilter 7

# Create a summary report
GALEON_Report.py -clust 3_OneFam_PhysEvoDistances_GFF3/ -plots Plots -ssize ChrSizes.txt -echo False
```

### 4.4. Test 4. Single family analysis using physical and evolutionary distances (using a protein MSA)
- Input directory with annotation files: `-a GFFs`
- Input g value: `-g 100,200`
- Enable the usage of Proteins: `-e enabled`
- Protein directory with fasta files: `-p Proteins`
- Use pre-computed MSA files: `-pm True`
- Output directory: `-outdir 4_OneFam_PhysEvoDistances_GFF3_pm`
- Chromosome/Scaffold size file: `-ssize ChrSizes.txt`
- Summarize the results for the first 7 largest scaffolds: `-sfilter 7`
- Write the path to files in the final report?: `-echo False`

```
cd Test_4

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
```

### 4.5. Test 5. Joint analysis of two gene families
- Input directory with annotation files: `-a GFFs`
- Input g value: `-g 100`
- Enable the usage of Proteins: `-e enabled`
- Color map: `-cmap_1 Blues_r`
- Square frame color for family 1 `-f red`
- Square frame color for family 2 `-f2 orange`
- Perform the joint analysis of two families: `-F BetweenFamilies`
- Output directory: `-outdir 5_TwoFamJointAnalysis_GFF3`

```
cd Test_5

# Run Galeon
GALEON_ControlScript.py clusterfinder -a GFFs/ -F BetweenFamilies -e disabled -outdir 5_TwoFamJointAnalysis_GFF3 -cmap_1 Blues_r -f red -f2 orange

# Generate summary files and tables
GALEON_SummaryFiles.py -fam merged -clust 5_TwoFamJointAnalysis_GFF3/ -coords GFFs/merged_dir/ -ssize ChrSizes.txt -sfilter 7

# Create a summary report
GALEON_Report.py -clust 5_TwoFamJointAnalysis_GFF3/ -ssize ChrSizes.txt -echo False
```

### 4.6. Test 6. GR and IR, single family analysis using physical and evolutionary distances (using a protein MSA) in *Dysdera silvatica*

This dataset includes the information for 98 GRs and 411 IRs annotated in *Dysdera silvatica* genome ([Escuer et al. 2022](https://doi.org/10.1111/1755-0998.13471)).

- Input directory with annotation files: `-a GFFs`
- Input g value: `-g 100`
- Enable the usage of Proteins: `-e enabled`
- Protein directory with fasta files: `-p Proteins`
- Use pre-computed MSA files: `-pm True`
- Output directory: `-outdir 4_OneFam_PhysEvoDistances_GFF3_pm`
- Square frame color: `-f orange`
- Chromosome/Scaffold size file: `-ssize ChrSizes.txt`
- Summarize the results for the first 7 largest scaffolds: `-sfilter 7`
- Write the path to files in the final report?: `-echo False`

```
cd Test_6

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
```

## 5. Citation

Vadim Pisarenco, Joel Vizueta, Julio Rozas. GALEON: A Comprehensive Bioinformatic Tool to Analyse and Visualise Gene Clusters in Complete Genomes. Submitted. 2024. https://www.biorxiv.org/content/10.1101/2024.04.15.589673v1

## 6. References

[1] Escuer, P., Pisarenco, V. A., Fernández-Ruiz, A. A., Vizueta, J., Sánchez-Herrero, J. F., Arnedo, M. A., Sánchez-Gracia, A., & Rozas, J. (2022). The chromosome-scale assembly of the Canary Islands endemic spider Dysdera silvatica (Arachnida, Araneae) sheds light on the origin and genome structure of chemoreceptor gene families in chelicerates. Molecular Ecology Resources, 22, 375–390. https://doi.org/10.1111/1755-0998.13471 

## 7. Troubleshooting

Should you encounter any error, please create an issue on GitHub specifying the error and providing as many details as possible.


