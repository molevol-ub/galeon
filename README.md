# GALEON

#### A Comprehensive Bioinformatic Tool to Analyse and Visualise Gene Clusters in Complete Genomes


To facilitate the identification, analysis, and visualisation of physically clustered gene family genes within chromosome-level genomes, we introduce GALEON, a user-friendly bioinformatic tool. GALEON identifies gene clusters by studying the spatial distribution of pairwise physical distances among gene family members along with the genome-wide gene density. The pipeline also enables the simultaneous analysis and comparison of two gene families, and allows the exploration of the relationship between physical and evolutionary distances. This tool offers a novel approach for studying the origin and evolution of gene families.


GALEON documentation can be also be found in: http://www.ub.edu/softevol/galeon


### Version history

V1: Initial release

## 0. Contents

 1. Installation and prerequisites
 2. Running GALEON
 3. Output
 4. Example data
 5. Citation
 6. Troubleshooting


## 1. Installation and prerequisites

GALEON is distributed as a set of several scripts that can be called from Galeon_masterScripts folder, but do not require any specific installation or compiilation step. However, the pipelines does require the following python modules and R packages.

### 1.1 - Python Packages
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

### 1.2 - R Packages

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

### 1.3 - Additional software

The following programs must be installed and available from command line: **mafft**, **bedtools**, **FastTree** and **iqtree2**. 

We provide a `bin` directory with binaries of `bedtools`, `FastTree` and `iqtree2`. If the Galeon conda environment was created, `mafft` and `newick_utils` should be available upon environment activation. Alternatively, check the corresponding documentation for installation instructions.


**Notes**:
- **Binaries** - Remember to make them executable by running: `chmod +x bin/*`.
- **FastTree** - By default, GALEON uses FastTree to infer the phylogeny among gene family copies. http://www.microbesonline.org/fasttree/
- **Iqtree** - IqTree can be used instead of FastTree to reconstruct the gene family phylogeny: http://www.iqtree.org/



## 2. Running GALEON

### 2.1. Estimate g parameter

The following script will estimate the g parameter based in the inputs:

-n define ...
-s ...
-g ...

```
python ControlScript.py gestimate -n 411 -s 1365.69 -g 100
```



### 2.2 Identifying gene clusters




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


