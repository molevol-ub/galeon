import os, subprocess, sys, ast, argparse, shutil

temp = subprocess.run("which GALEON_SummaryFiles.py", shell=True, capture_output=True, text=True).stdout.strip()
PATHtoGaleonScripts = os.path.dirname(temp)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Summary Plots', formatter_class=argparse.RawTextHelpFormatter) #ArgumentDefaultsHelpFormatter)
    parser._optionals.title = "Input arguments"
    
    parser.add_argument("-fam", "--FamilyName",
        type = str,
        required = True,
        help = "(REQUIRED). Gene family name (i.e: 'GR', 'IR' or 'merged' when studying two families at the same time.")

    parser.add_argument("-clust", "--ClusterfinderResults",
        type = str,
        required = True,
        help = "(REQUIRED). GALEON clusterfinder results directory.")

    parser.add_argument("-coords", "--CoordsDir",
        type = str,
        required = True,
        help = "(REQUIRED). Coordinates directory of the input family.")

    parser.add_argument("-ssize", "--ScaffoldSizeFile",
        type = str,
        required = True,
        help = "(REQUIRED). Chromosome/Scaffold size file. A tab-separated file composed of three columns: 'ScaffoldID', 'Scaffold length', 'Scf. Associated name'.")
    
    parser.add_argument("-sfilter", "--ScaffoldToRetain",
        type = str,
        required = True,
        help = "(REQUIRED). Scaffolds to retain. Choose \n-'ALL' - retain all of them. \
                                                        \n 'NUM' - retain only the firts NUM largest scaffodls. \
                                                        \n 'Provide a FILE' - retain only some scaffolds of interest. File format: a list of scaffolds in one column")

    parser.add_argument("-colopt", "--ColorOption",
        type = str,
        default = "Palette",
        choices = ['OneColor', 'Palette'],
        help = "Color option for Summary Plots")

    parser.add_argument("-colval", "--ColorValue",
        type = str,
        default = "mako",
        help = "Color or color palette for Summary Plots.")

    parser.add_argument("-outdir", "--OutDir",
        type = str,
        default = "Plots",
        help = "Output directory to store the summary plots")    
    # Parsing arguments
    args = parser.parse_args()


    # Get configvalue
    config = vars(args)
    print(config)


# In[18]:


# Input results directory
# RESULTSdir = "clusterfinder_Results_Directory/"
# COORDSdir = "GFFs/"
# CHRSizefile = "ChrSizes.txt"
# GENEFAMILY = "IR"
# SCFtoretain = "3" # ALL; NUM; File. If a Num, provide as string

# RESULTSdir = "clusterfinder_Results_Directory_mm/"
# COORDSdir = "GFFs/merged_dir/"
# CHRSizefile = "ChrSizes.txt"
# GENEFAMILY = "merged"
# SCFtoretain = "7" # ALL; NUM; File. If a Num, provide as string

# Input
RESULTSdir = config["ClusterfinderResults"]
COORDSdir = config["CoordsDir"]
CHRSizefile = config["ScaffoldSizeFile"]
GENEFAMILY = config["FamilyName"]
SCFtoretain = config["ScaffoldToRetain"]
ColorOpt = config["ColorOption"]
ColorValue = config["ColorValue"]
PlotsOutDir =  config["OutDir"]


if os.path.exists(f"{RESULTSdir}/{PlotsOutDir}"):
    pass
else:
    os.mkdir(f"{RESULTSdir}/{PlotsOutDir}")

# Don't modify the below code
# Script to create the output table
ScriptsDir = f"{PATHtoGaleonScripts}/Scripts"
S_script = f"{ScriptsDir}/CreateDataframe.py"

if GENEFAMILY == "merged":
    PlotScript = f"{ScriptsDir}/ClusterSize_Figure_merged.py"
else:
    PlotScript = f"{ScriptsDir}/ClusterSize_Figure.py"

# RESULTSdir = sys.argv[1]
# COORDSdir = sys.argv[2]
# CHRSizefile = sys.argv[3]
# GENEFAMILY = sys.argv[4]


# In[19]:


def load_dict(ifile):
    with open(ifile) as f1:
        temp = f1.read()

        D_clusters = ast.literal_eval(temp)
    return D_clusters

def CheckClusterData(i_file):
    # Load dictionary
    D_clust = load_dict(i_file)
    
    # Save here the size of clusters
    clustinfo_count = set()

    for k,v in D_clust.items():
        clustinfo_count.add(len(v))
    
    # If there are no clusters...
    if clustinfo_count == {0}:
        print(f"No clusters data in this dictionary: '{i_file}'")
        return None
    else:
        return D_clust


# In[20]:


PlotsOutDir = f"{GENEFAMILY}_fam"

if os.path.exists(PlotsOutDir) == False:
    
    os.mkdir(PlotsOutDir)
else:
    shutil.rmtree(PlotsOutDir)
    os.mkdir(PlotsOutDir)
    
# In[21]:


# Find cluster dictionaries
ClusterDictList = []
for r, d, f in os.walk(RESULTSdir):
    if "_matrices_" in r and f"{GENEFAMILY}_" in os.path.basename(r):
        
        temp = [i for i in f if i.endswith(".cluster.dict")]
        
        # checkpoint
        if len(temp) == 0:
            raise ValueError("No cluster dictionary file was detected")
        elif len(temp) > 1:
            raise ValueError("More than one dictionary file was detected. Only one should be present in the input directory")
        elif len(temp) == 1:
            DictFile = temp[0]
            DictFilePATH = f"{r}/{DictFile}"
            ClusterDictList.append(DictFilePATH)
        else:
            raise ValueError("Unknown error")
            
# checkpoint
if len(ClusterDictList) == 0:
    emsg = "No Cluster dictionary was found"
    raise ValueError(emsg)


# In[22]:


''' Find the input coordinates file '''

CoordsFileList = [i for i in os.listdir(COORDSdir) if i.endswith("collapsed.temp") and i.startswith(f"{GENEFAMILY}_")]
if len(CoordsFileList) > 1:
    emsg = f"More than one coordinate file has been found for this family '{GENEFAMILY}'"
    raise ValueError(emsg)
elif len(CoordsFileList) == 0:
    emsg = "No coordinates file was found!"
    raise ValueError(emsg)
else:
    CoordsFile = CoordsFileList[0]

''' Find the input coordinates file format '''
CoordsFileFormat = ""

for i in [".gff3", ".bed1", ".bed2"]:
    if i in CoordsFile:
        CoordsFileFormat = i.replace(".","")
        
# checkpoint
if CoordsFileFormat == None:
    raise ValueError("Input Coordinate File has an unknown format")
else:
    pass

CoordsFile = f"{COORDSdir}/{CoordsFile}"


# In[23]:


''' Check for the presence of ChrSizes file '''
if os.path.exists(CHRSizefile):
    pass
else:
    emsg = f"Chromosome/Scaffold sizes file '{CHRSizefile}' doesn't exist"
    raise ValueError(emsg)


# In[24]:


''' Run Create DataFrame script '''

for iClustDictFile in ClusterDictList:
    if CheckClusterData(iClustDictFile) == None:
        pass
    else:
        gvalue_hint = os.path.basename(iClustDictFile).split("_matrices_")[1].split(".Dict")[0]
        output_table = f"{GENEFAMILY}_family_GeneLocation.table.{gvalue_hint}.tsv"
        # print(CoordsFile, iClustDictFile, CHRSizefile, output_table)

        # print(" ".join(["python", S_script, CoordsFile, iClustDictFile, CHRSizefile, f"{PlotsOutDir}/{output_table}", SCFtoretain]))
        temp = subprocess.run(["python", S_script, CoordsFile, iClustDictFile, CHRSizefile, f"{PlotsOutDir}/{output_table}", SCFtoretain], text=True, capture_output=True, check=True)
        temp = subprocess.run(["python", PlotScript, f"{PlotsOutDir}/{output_table}", gvalue_hint, PlotsOutDir, SCFtoretain, ColorOpt, ColorValue], text=True, capture_output=True, check=True)


shutil.move(PlotsOutDir, f'{config["ClusterfinderResults"]}/{config["OutDir"]}')
