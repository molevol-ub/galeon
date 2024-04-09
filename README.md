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

GALEON is distributed as a set of scripts that can be called from Galeon_masterScripts folder, but do not require any specific installation or compilation step. However, the pipeline does require the following python modules, R packages, as well as some external software.

### 1.1. Python Packages
If conda is available, all these packages can be easily obtained by installing the Galeon conda environment (see 1.4), as well as some of the required software : mafft and newick_utils. Alternatively, you may install them separately using `pip`, consult the appropriate documentation for each of them.

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
conda env create -f GaleonEnv.yml

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
- 4 tab separated columns
- This provides directly the coordinates and the gene IDs.

| Scaffold ID | start | end | attribute |
| :-------------: | :-------------: | :-------------: | :-------------: |
| Scaffold_14804_HRSCAF_18385 | 41841903 | 41843055 | g10232 |
| Scaffold_14804_HRSCAF_18385 | 47268322 | 47268742 | g10331 |
| Scaffold_14804_HRSCAF_18385 | 47277448 | 47277868 | g10332 |


#### 2.1.3. BED1 format
- 3 tab separated columns
- This format does not include gene IDs. Gene names will be given by the gene order. Consider this example annotation file name: `GR_fam.bed1`. The genes will be named as *GR_0*, *GR_1*, etc.


| Scaffold ID | start | end |
| :-------------: | :-------------: | :-------------: |
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
| :-------------: | :-------------: | :-------------: |
| Scaffold_15362_HRSCAF_19823 | 317950935 | ChrX |
| Scaffold_14804_HRSCAF_18385 | 177171321 | Chr1 |
| Scaffold_14178_HRSCAF_16784 | 176727214 | Chr2 |

<br>

## 3. Running GALEON

### 3.1. Estimate g parameter (*mode: gestimate*)
In this mode, the pipeline estimates the expected number of genes found in a number of bases, as well as the number of genes expected across the g input values and the probability of finding 2 or more genes in a window of g size (i.e.: 100 Kb), which would be considered as a cluster in the following analyses (**Section 3.2.**).

Run the following command to estimate the g parameter based on the inputs. No input files are required here.

- `-n NUM` Gene family size (number of genes)
- `-s NUM` Genome size (in Mb units)
- `-g NUM` g value (in Kb units): more that one can be tested
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
2. column: **Exp. 1 gene each X Mb** - 1 gene, expected to be found each "X" Mb
3. column: **Exp. genes / Mb** - # of genes, expected to be found each Mb
4. column: **Exp. Genes / Kb** - # of genes, expected to be found each 1 Kb
5. column: **Exp. Genes/g value** - # of genes, expected to be found each "g" Kb
6. column: **P(X>=2) / g value** - The probability of finding by chance two (or more) genes in a "g" kb stretch
7. column: **Poisson's λ**

<br>

### 3.2. Gene cluster identification (*mode: clusterfinder*)
In this mode the pipeline analyzes one (or several) gene families to identify clusters of genes in the genome. 

**Help message**
```
GALEON_ControlScript.py gestimate -h
```

#### 3.2.1. One family analysis using physical distances
**Inputs**

- More than one families can be analyzed at the same time

```
# How your input annotation directory "GFFs/" should look
└── GFFs
    ├── GR_fam.gff3
    └── IR_fam.gff3
```

**Commands**

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
    - This directory contains summary plots with cluster size distribution at genome level (dir.: `SummaryPlots_1.0g`) and individual scaffold level (dir.: `IndividualPlots_1.0g`)
    - Also, it contains some summary tables.
        - `X_family_ClusterSizes.table.1.0g.tsv` - This tables contains a more explicit information about the size of each cluster in each scaffold.
        - `X_family_GeneLocation.table.1.0g.tsv` - This tables contains information detailed at gene level, specifically, gene coordinates and membership (singleton or clustered).
        - `X_family_GeneOrganizationGenomeSummary.table.1.0g.tsv` - This tables informs about how many genes are organized in clusters (Clustered category) or are not (Singleton category) at genome level.
        - `X_family_GeneOrganizationSummary.table.1.0g.tsv` - This tables informs about how many genes are organized in clusters (Clustered category) or are not (Singleton category) in each scaffold.

- **Reports** directory
    - This directory contains html reports, one for each family and g value.


Galeon results directory, example tree-like representation for two gene families "GR_fam" and "IR_fam":

```
clusterfinder_Results_Directory/
├── PhysicalDist_Matrices
│   ├── GR_fam.gff3.temp_matrices # Physical distance matrices, *matrix (in bp units)
│   ├── IR_fam.gff3.temp_matrices 
│   ├── GR_fam.gff3.temp_matrices_1.0g # Heatmaps in svg and pdf format
│   └── IR_fam.gff3.temp_matrices_1.0g
|
├── Plots # Contains summary tables and plots for each input family
│   ├── GR_fam
│   │   ├── GR_family_ClusterSizes.table.1.0g.tsv
│   │   ├── GR_family_GeneLocation.table.1.0g.tsv
│   │   ├── GR_family_GeneOrganizationGenomeSummary.table.1.0g.tsv
│   │   ├── GR_family_GeneOrganizationSummary.table.1.0g.tsv
│   │   |
│   │   ├── IndividualPlots_1.0g 
│   │   └── SummaryPlots_1.0g
│   └── IR_fam
│   │   ├── IR_family_ClusterSizes.table.1.0g.tsv
│   │   ├── IR_family_GeneLocation.table.1.0g.tsv
│   │   ├── IR_family_GeneOrganizationGenomeSummary.table.1.0g.tsv
│   │   ├── IR_family_GeneOrganizationSummary.table.1.0g.tsv
│   │   |
│   │   ├── IndividualPlots_1.0g
│   └── └── SummaryPlots_1.0g
|
└── Reports # One report for each family and tested g value.
    ├── GR_fam_1.0g_Report.html
    └── IR_fam_1.0g_Report.html
```

<br>

#### 3.2.2. One family analysis using physical and evolutionary distances
**Inputs**

- More than one families can be analyzed at the same time

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

**Step 2-3)** Generate summary plots, tables and the HTML report

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
    - This directory containts pairwise physical + evolutionary distance files (`*matrix` and `*distance.tsv`) between genes of a given Gene Family, in each scaffold.
    - Also, it contains physical + evolutionary distance heatmaps and scatterplots in svg and pdf format at a genome level (representing the entire gene family) as well as for each scaffold (showing each one separately).
- **PhysicalDist_Matrices** directory. Same as in Section 3.2.1.
- **Plots** directory. Same as in Section 3.2.1.
- **Reports** directory. Same as in Section 3.2.1.

Galeon results directory, example tree-like representation for two gene families "GR_fam" and "IR_fam":

```
clusterfinder_Results_Directory/
├── GR_fam.GlobScatterPlot_1.0g.pdf # Physical vs Evo. distance scatter plot considering at genome level, that is, considering all the genes of the input gene family
├── GR_fam.GlobScatterPlot_1.0g.svg
├── IR_fam.GlobScatterPlot_1.0g.pdf
├── IR_fam.GlobScatterPlot_1.0g.svg
│
├── MergedDistances_Dataframes
│   ├── GR_fam.IntermediateFiles
│   ├── GR_fam.merged.matrices # Physical + Evolutionary distance matrices
│   ├── GR_fam.plots_1.0g # Physical + Evolutionary distance Heatmaps and Scatterplots in svg and pdf format
│   ├── IR_fam.IntermediateFiles
│   ├── IR_fam.merged.matrices
│   └── IR_fam.plots_1.0g
│           
├── PhysicalDist_Matrices 
│   ├── GR_fam.gff3.temp_matrices  # Physical distance matrices, *matrix (in bp units)
│   ├── GR_fam.gff3.temp_matrices_1.0g # Heatmaps in svg and pdf format
│   ├── IR_fam.gff3.temp_matrices
│   └── IR_fam.gff3.temp_matrices_1.0g
│
├── Plots # same content as in "Section 3.2.1"
│   ├── GR_fam
│   │   ├── IndividualPlots_1.0g
│   │   └── SummaryPlots_1.0g
│   └── IR_fam
│       ├── IndividualPlots_1.0g
│       └── SummaryPlots_1.0g
│
└── Reports # same content as in "Section 3.2.1"
    ├── GR_fam_1.0g_Report.html
    └── IR_fam_1.0g_Report.html

```

#### 3.2.3. Joint family analysis
Find clusters between two input families using the coordinates from the input files. ONLY two families can be analysed at once, note that Proteins cannot be used in this mode.

**Inputs**

- Only two families are allowed and can be analyzed at the same time.

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

Note how the `-coord` parameter is specified here (`-coord GFFs/merged_dir`), it is a bit different because you must add "merged_dir" at the end of the command.

Galeon results directory content

- **PhysicalDist_Matrices** directory
    - This directory containts pairwise physical distance matrices (`*matrix`) between genes of a given Gene Family, in each scaffold.
    - Also, physical distance heatmaps are present in svg and pdf format.
- **Plots** directory
    - This directory contains summary plots with cluster size distribution at genome level (dir.: `SummaryPlots_1.0g`) and individual scaffold level (dir.: `IndividualPlots_1.0g`)
    - Since now we are interested in the joint analysis of two families (let's name them X and Y), the output summary plots are presented in 4 subdirectories:
        - `X` - plots for clusters of X gene family
        - `Y` - plots for clusters of Y gene family
        - `X.Y` - plots for clusters of X and Y gene family members (when members both families belong to the same cluster).
        - `merged` - plots for clusters of X and Y gene family considered as one.
    - Also, it contains some summary tables.
        - `X_family_ClusterSizes.table.1.0g.tsv` - This tables contains a more explicit information about the size of each cluster in each scaffold.
        - `merged_family_GeneLocation.table.1.0g.tsv` - This tables contains information detailed at gene level of both families, specifically, gene coordinates and membership (singleton or clustered).
        - `X_family_GeneOrganizationGenomeSummary.table.1.0g.tsv` - This tables informs about how many genes are organized in clusters (Clustered category) or are not (Singleton category) at genome level.
        - `X_family_GeneOrganizationSummary.table.1.0g.tsv` - This tables informs about how many genes are organized in clusters (Clustered category) or are not (Singleton category) in each scaffold.
        - `X.Y_family_GeneOrganizationSummary.table.1.0g.tsv` - This tables informs about how many genes of both families are jointly organized in clusters (Clustered category) or are not (Singleton category) in each scaffold. This file appears only when there are clusters formed by members of both families. Otherwise, it is absent.

- **Reports** directory
    - This directory contains html reports, one for each family and g value.

Galeon results directory, example tree-like representation for two gene families "GR_fam" and "IR_fam":

```
clusterfinder_Results_Directory/
├── PhysicalDist_Matrices
│   ├── merged_fam.gff3.temp_matrices # contains matrices *matrix
│   └── merged_fam.gff3.temp_matrices_1.0g # physical distance plots in svg and pdf format
│
├── Plots
│   └── merged_fam # summary tables and plots
│       ├── GR_family_ClusterSizes.table.1.0g.tsv
│       ├── GR_family_GeneOrganizationGenomeSummary.table.1.0g.tsv
│       ├── GR_family_GeneOrganizationSummary.table.1.0g.tsv
│       │
│       ├── IR_family_ClusterSizes.table.1.0g.tsv
│       ├── IR_family_GeneOrganizationGenomeSummary.table.1.0g.tsv
│       ├── IR_family_GeneOrganizationSummary.table.1.0g.tsv
│       │
│       ├── GR.IR_family_GeneOrganizationSummary.table.1.0g.tsv
│       │
│       ├── merged_family_GeneLocation.table.1.0g.tsv
│       │
│       ├── IndividualPlots_1.0g
│       │   ├── GR # Gr family clusters' size distribution
│       │   ├── GR.IR # Two family clusters size distribution
│       │   ├── IR # Gr family clusters' size distribution
│       │   └── merged # All clusters' size distribution
│       └── SummaryPlots_1.0g
│           ├── GR
│           ├── GR.IR
│           ├── IR
│           └── merged
│           
└── Reports
    ├── merged_fam_1.0g_Report.html
    └── merged_fam_1.0g_Report.Rmd
```

<br>


## 4. Example dataset

Severals tests are available in [Example_data](https://github.com/molevol-ub/galeon/tree/main/Example_data) directory. You enter to each *test* directory and run the commands described below. All the example are solved, check this file [Solved_examples.tar.gz](https://github.com/molevol-ub/galeon/blob/main/GALEON_compressed_file/Solved_examples.tar.gz)
- Remember to have the Galeon conda environment activated before running the commands.


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

### 4.2. Test 2. One family analysis using physical distances

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

### 4.3. Test 3. One family analysis using physical and evolutionary distances (use Proteins)
- Input directory with annotation files: `-a GFFs`
- Input g value: `-g 100,200`
- Enable the usage of Proteins: `-e enabled`
- Protein directory with fasta files: `-p Proteins`
- Output directory: `-outdir 3_OneFam_PhysEvoDistances_GFF3`
- Chromosome/Scaffold size file: `-ssize ChrSizes.txt`
- Summarize the results for the first 7 largest scaffolds: `-sfilter 7`
- Write the path to files in the final report?: `-echo False`

```
cd Test_3

# Run Galeon
GALEON_ControlScript.py clusterfinder -a GFFs/ -g 100,300 -e enabled -p Proteins/ -outdir 3_OneFam_PhysEvoDistances_GFF3/

# Get evolutionary statistics (Cst) and perform the Mann-Whitney test
GALEON_GetEvoStats.py -clust 3_OneFam_PhysEvoDistances_GFF3/ -prot Proteins/ -coords GFFs

# Generate summary files and tables for GR family
GALEON_SummaryFiles.py -fam GR -clust 3_OneFam_PhysEvoDistances_GFF3/ -coords GFFs -ssize ChrSizes.txt -sfilter 7

# Generate summary files and tables for IR family
GALEON_SummaryFiles.py -fam IR -clust 3_OneFam_PhysEvoDistances_GFF3/ -coords GFFs -ssize ChrSizes.txt -sfilter 7

# Create a summary report
GALEON_Report.py -clust 3_OneFam_PhysEvoDistances_GFF3/ -plots Plots -ssize ChrSizes.txt -echo False
```

### 4.4. Test 4. One family analysis using physical and evolutionary distances (use MSA)
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
- Input g value: `-g 100,200`
- Enable the usage of Proteins: `-e enabled`
- Perform the joint analysis of two families: `-F BetweenFamilies`
- Output directory: `-outdir 5_TwoFamJointAnalysis_GFF3`

```
cd Test_5

# Run Galeon
GALEON_ControlScript.py clusterfinder -a GFFs/ -F BetweenFamilies -e disabled -outdir 5_TwoFamJointAnalysis_GFF3

# Generate summary files and tables
GALEON_SummaryFiles.py -fam merged -clust 5_TwoFamJointAnalysis_GFF3/ -coords GFFs/merged_dir/ -ssize ChrSizes.txt -sfilter 7

# Create a summary report
GALEON_Report.py -clust 5_TwoFamJointAnalysis_GFF3/ -ssize ChrSizes.txt -echo False
```


## 5. Citation

Vadim Pisarenco, Joel Vizueta, Julio Rozas. GALEON . Submitted. 2023.


## 6. Troubleshooting

Should you encounter any error, please create an issue on GitHub specifying the error and providing as many details as possible.


