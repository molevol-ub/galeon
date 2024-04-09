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


## 2. Running GALEON

### 2.1. Estimate g parameter

The following script will estimate the g parameter based in the inputs:

-n define ...
-s ...
-g ...

```
python ControlScript.py gestimate -n 411 -s 1365.69 -g 100
```



### 2.2. Identifying gene clusters




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


