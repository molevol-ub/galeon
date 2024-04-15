#!/home/vadim/miniconda3/envs/Galeon/bin/python
import os, re, subprocess, shutil, argparse
import pandas as pd

temp = subprocess.run("which GALEON_Report.py", shell=True, capture_output=True, text=True).stdout.strip()
PATHtoGaleonScripts = os.path.dirname(temp)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Final Report', formatter_class=argparse.RawTextHelpFormatter) #ArgumentDefaultsHelpFormatter)
    parser._optionals.title = "Input arguments"
    
    parser.add_argument("-clust", "--ClusterfinderResults",
        type = str,
        required = True,
        help = "(REQUIRED). GALEON clusterfinder results directory.")
    
    parser.add_argument("-plots", "--PlotsDir",
        type = str,
        default = "Plots",
        help = "Input plots directory.")

    parser.add_argument("-MW", "--MWDir",
        type = str,
        default = "MannWhitney_StatisticsResults",
        help = "Mann-Whitney results directory")

    parser.add_argument("-ssize", "--ScaffoldSizeFile",
        type = str,
        required = True,
        help = "(REQUIRED). Chromosome/Scaffold size file. A tab-separated file composed of three columns: 'ScaffoldID', 'Scaffold length', 'Scf. Associated name'.")
    
    parser.add_argument("-echo", "--EchoPaths",
        type = str,
        default = "False",
        choices = ["True", "False"],
        help = "Print paths to the graphs in the final report")

    parser.add_argument("-outdir", "--ReportDir",
        type = str,
        default = "Reports",
        help = "Output report directory (Default: Reports)")
        
    # Parsing arguments
    args = parser.parse_args()


    # Get configvalue
    config = vars(args)
    print(config)


# In[ ]:


# Input dirs
GaleonClusterFinderResults = config["ClusterfinderResults"]
PLOTSdir = config["PlotsDir"]
PLOTSdir = f"{GaleonClusterFinderResults}/{PLOTSdir}"
MWtest = config["MWDir"]
MWtest = f"{GaleonClusterFinderResults}/{MWtest}"

# Scaffold sizes file
CHRsizeFile = config["ScaffoldSizeFile"]

# Outdir
OUTdir = config["ReportDir"]

echo_paths =config["EchoPaths"]

if echo_paths == "True":
    echo_paths = True
elif echo_paths == "False":
    echo_paths = False


# Mann-Whitney data?
if os.path.exists(PLOTSdir) == False:
    emsg = f"ERROR! This directory doesn't exist: '{PLOTSdir}'"
    raise ValueError(emsg)
    
if os.path.exists(MWtest) == False:
    emsg = f"WARNING! This directory doesn't exist: '{MWtest}'"
    MannWhitney_results = False
    EvoDistUse = False
else:
    MannWhitney_results = True
    # Use of evolutionary distance
    EvoDistUse = True


# In[ ]:


'''Guide list of chr/scf sizes'''

# Load chr size file
DF_ChrSize = pd.read_csv(CHRsizeFile, sep="\t", header=None)
DF_ChrSize.columns = ["ScaffoldID", "Scaffold length", "ScfAssociatedName"]
DF_ChrSize = DF_ChrSize.sort_values(by="Scaffold length", ascending=False)

# Create guidelists
guidelist = [] # Create a guidelist to sort the output individual plots by scf associated names
scfguidelist = DF_ChrSize["ScaffoldID"].to_list() # Create a guidelist to sort the output individual plots by scaffold IDs

for i in DF_ChrSize["ScfAssociatedName"].to_list(): 
    if i not in guidelist:
        guidelist.append(i)


# Create the output directory
if os.path.exists(OUTdir) == False:
    os.mkdir(OUTdir)
else:
    shutil.rmtree(OUTdir)
    os.mkdir(OUTdir)
    pass


# Dictionary that will store the paths to the figures and summary files
D_GLOB = {}

# Find summary tables
SummaryTable_list = []
for r, d, f in os.walk(PLOTSdir):
    if len(f) != 0:
        for ifile in f:
            if "_family_" in ifile and ifile.endswith(".tsv"):
                path2table = f"{r}/{ifile}"
                SummaryTable_list.append(path2table)

''' Find summary files '''
for sumfile in SummaryTable_list:
    GeneFamName = os.path.dirname(sumfile).split("/")[-1]
    gvalue_key = sumfile.split(".table.")[-1].replace(".tsv","")
    # print(GeneFamName, gvalue_key, sumfile)
    
    if GeneFamName not in D_GLOB.keys():
        D_GLOB[GeneFamName] = {}
        D_GLOB[GeneFamName][gvalue_key] = {"SummaryFiles" : [sumfile]}
    else:
        if gvalue_key not in D_GLOB[GeneFamName].keys():
            D_GLOB[GeneFamName][gvalue_key] = {"SummaryFiles" : [sumfile]}
        else:
            D_GLOB[GeneFamName][gvalue_key]["SummaryFiles"].append(sumfile)

''' Find Physical distance plots '''
for r, d, f in os.walk(GaleonClusterFinderResults):
    if len(f) != 0:
        # Physical heatmaps
        filelist = [f"{r}/{i}" for i in f if "PhysicalDist_Heatmap" in i and i.endswith(".svg")]
        if len(filelist) != 0:
            FamName = r.split("/")[-1].split("_fam")[0] + "_fam"
            gvalue = r.split("_")[-1]
            
            try:
                D_GLOB[FamName][gvalue]["GaleonClusterFinderResults"] = filelist
            except KeyError:
                pass
        
        if EvoDistUse == True:
            # Physical + Evo heatmaps
            # filelist = [f"{r}/{i}" for i in f if "PhysicalEvolutionaryDist" in i and i.endswith(".svg")]
            filelist = [i for i in f if "PhysicalEvolutionaryDist" in i and i.endswith(".svg")]
            if len(filelist) != 0:
                FamName = r.split("/")[-1].split("_fam")[0] + "_fam"
                gvalue = r.split("_")[-1]
                
                templist = ["_".join(i.split("_")[:-2]) for i in filelist]
                minidict = {i : [f"{r}/{i}_Heatmap_{gvalue}.svg", f"{r}/{i}_Scatterplot_{gvalue}.svg"] for i in templist}
                
                try:
                    D_GLOB[FamName][gvalue]["PhysicalEvoDistance_Figures"] = minidict
                except KeyError:
                    pass


''' Find Summary plots '''
for r, d, f in os.walk(PLOTSdir):

    if len(d) != 0 and len(f) != 0:
        
        IndividualPlots_dirs = [i for i in d if "IndividualPlots" in i]
        SummaryPlots_dirs = [i for i in d if "SummaryPlots" in i]
        
        # Add IndividualPlots
        if len(IndividualPlots_dirs) != 0:
            for idir in IndividualPlots_dirs:
                temp = f"{r}/{idir}"
                GeneFamname = os.path.basename(r)
                gvalue_Key = idir.split("_")[-1]
                
                path2files = [f"{temp}/{f}" for f in os.listdir(temp)]
                D_GLOB[GeneFamname][gvalue_Key][idir] = path2files
        
        # Add SummaryPlots
        if len(SummaryPlots_dirs) != 0:
            for idir in SummaryPlots_dirs:
                temp = f"{r}/{idir}"
                
                GeneFamname = os.path.basename(r)
                gvalue_Key = idir.split("_")[-1]
                
                path2files = [f"{temp}/{f}" for f in os.listdir(temp)]
                D_GLOB[GeneFamname][gvalue_Key][idir] = path2files
                
if EvoDistUse == True:     
    for r, d, f in os.walk(GaleonClusterFinderResults):
        if len(f) != 0:

            # Physical + Evo heatmaps
            # filelist = [f"{r}/{i}" for i in f if "PhysicalEvolutionaryDist" in i and i.endswith(".svg")]
            filelist = [i for i in f if "PhysicalEvolutionaryDist" in i and i.endswith(".svg")]
            if len(filelist) != 0:
                FamName = r.split("/")[-1].split("_fam")[0] + "_fam"
                gvalue = r.split("_")[-1]

                templist = ["_".join(i.split("_")[:-2]) for i in filelist]
                minidict = {i : [f"{r}/{i}_Heatmap_{gvalue}.svg", f"{r}/{i}_Scatterplot_{gvalue}.svg"] for i in templist}

                try:
                    D_GLOB[FamName][gvalue]["PhysicalEvoDistance_Figures"] = minidict
                except KeyError:
                    pass

            globplots_filelist = [i for i in f if "GlobScatterPlot" in i and i.endswith(".svg")]
            if len(globplots_filelist) != 0:
                for ifile in globplots_filelist:
                    FamName = ifile.split("/")[-1].split("_fam")[0] + "_fam"
                    gvalue = ifile.split("_")[-1].replace(".svg","")
                    sumkey = f"SummaryPlots_{gvalue}"
                    D_GLOB[FamName][gvalue][sumkey].append(f"{r}/{ifile}")

# In[ ]:


''' Find Mann-Whitney results '''
if MannWhitney_results == True:
    for i in os.listdir(MWtest):
        temp = i.split("_")
        GeneFamID = "_".join(temp[:2])

        gValueKey = temp[-1].split(".")[-3:-1]
        gValueKey = ".".join(gValueKey)

        if "rawdata" in i:
            D_GLOB[GeneFamID][gValueKey]["MW raw data"] = f"{MWtest}/{i}"
        elif "results.brief" in i:
            D_GLOB[GeneFamID][gValueKey]["MW results"] = f"{MWtest}/{i}"
        elif "GlobalStats_value" in i:
            D_GLOB[GeneFamID][gValueKey]["Cst glob stats"] = f"{MWtest}/{i}"

# In[ ]:


# Fx: format scaffold names (before exporting to the report)
def AddChrNameLength(i_scaffold, i_guideDF):
    # Fx: add the associated chr/scf name and length
    
    temp = i_guideDF[i_guideDF["ScaffoldID"] == i_scaffold]
    if len(temp.index) != 1:
        raise ValueError("Unknown error")

    else:
        temp2 = temp.iloc[0].to_list()[::-1]
        temp2 = [str(i) for i in temp2]
        outname = " - ".join(temp2)

        return outname
    
def FigName_format(iname, iDF_ChrSize):
    # Fx: format scaffold name, the scf name is inferred directly from the filename
    
    temp = iname.split(".temp_matrices.")
    FamID = temp[0].split(".")[0]
    ScfID = temp[1].split(".matrix.")[0]
    
    ScfID = AddChrNameLength(ScfID, iDF_ChrSize)
    
    # outname = f"{FamID}: {ScfID}"
    outname = ScfID
    return outname



# In[ ]:


# Fx: writing functions to export the information

def write_table(i_path, o_file, i_opt, i_echo):
    if i_echo == False:
        i_echo = "FALSE"
    elif i_echo == True:
        i_echo = "TRUE"
        
    if i_opt == 1:
        text_list = ['```{' + 'r load-csv, echo={}'.format(i_echo) + '}','library(DT)',
        f'table_data <- read.csv("{i_path}", sep="\t")',
        'datatable(table_data, options = list(searching = TRUE, pageLength = 5))',
        '```\n']
    elif i_opt == 2:
        text_list = ['```{' + 'r, echo={}'.format(i_echo) + '}','library(DT)',
        f'table_data <- read.csv("{i_path}", sep="\t")',
        'datatable(table_data, options = list(searching = TRUE, pageLength = 5))',
        '```\n']
    
    for i in text_list:
        print(i, file=o_file, end="\n")

def write_PhysMxHeatmap_figures(i_dict, i_key, o_file, i_echo, i_opt=True):
    if i_echo == False:
        i_echo = "FALSE"
    elif i_echo == True:
        i_echo = "TRUE"
        
    filelist = Sort_IndividualPlots(i_dict[i_key], scfguidelist, "Singlefamily", True)
    for idx, fig in enumerate(filelist):
        idx += 1
        temp = os.path.basename(fig)
        figname = FigName_format(temp, DF_ChrSize)
        
        info = [f'**Figure {idx}. {figname}**', '<div>\n',
                '```{' + 'r, '+ 'out.width="75%", ' + 'echo={}'.format(i_echo) + '}',
                'library(knitr)',
                f'knitr::include_graphics("{fig}")',
                '```\n\n']
        for i in info:
            print(i, file=o_file, end="\n")
            
def write_figures(i_dict, i_key, o_file, i_echo, i_opt):
    if i_echo == False:
        i_echo = "FALSE"
    elif i_echo == True:
        i_echo = "TRUE"
        
    temp = i_dict[i_key]
    
    if i_opt == "SummaryPlots":
        figlist = sorted(temp, key=custom_sorter)
            
    elif i_opt == "IndividualPlots":
        figlist = Sort_IndividualPlots(temp, guidelist, "Singlefamily")
    
    for idx, fig in enumerate(figlist):
        idx += 1
        if "Global_ClusterSize" in fig:
            info = [f'**Figure {idx}. {os.path.basename(fig)}**', '<div>\n',
                    '```{' + 'r, ' + 'out.width="75%", '+ 'echo={}'.format(i_echo) + '}',
                    'library(knitr)',
                    f'knitr::include_graphics("{fig}")',
                    '```\n\n']
        elif "ClusterSize_Distribution" in fig:
            info = [f'**Figure {idx}. {os.path.basename(fig)}**', '<div>\n',
                    '```{' + 'r, ' + 'out.width="75%", '+ 'echo={}'.format(i_echo) + '}',
                    'library(knitr)',
                    f'knitr::include_graphics("{fig}")',
                    '```\n\n']
        
        else:
            info = [f'**Figure {idx}. {os.path.basename(fig)}**', '<div>\n',
                    '```{' + 'r, ' + 'out.width="65%", '+ 'echo={}'.format(i_echo) + '}',
                    'library(knitr)',
                    f'knitr::include_graphics("{fig}")',
                    '```\n\n']
        for i in info:
            print(i, file=o_file, end="\n")
            
def write_2_figures(i_dict, i_fam, i_gvalue, i_echo, o_file):
    if i_echo == False:
        i_echo = "FALSE"
    elif i_echo == True:
        i_echo = "TRUE"
        
    for idx, (kcase, vfigures) in enumerate(i_dict.items()):
        
        figname = FigName_format(kcase, DF_ChrSize)
        
        if idx == 0:
            fig1, fig2 = [*vfigures]

            info = [f'**Figure {idx+1}. {figname}**', '<div>\n',
                    '```{' + 'r image_grobs, fig.show = "hold", out.width = "50%", fig.align = "default", echo={}'.format(i_echo) + '}',
                    f'knitr::include_graphics("{fig1}")',
                    f'knitr::include_graphics("{fig2}")',
                    '```']
            for i in info:
                print(i, file=o_file, end="\n")
                
        else:
            
            fig1, fig2 = [*vfigures]

            info = [f'**Figure {idx+1}. {figname}**', '<div>\n',
                    '```{' + 'r fig.show = "hold", out.width = "50%", fig.align = "default", echo={}'.format(i_echo) + '}',
                    f'knitr::include_graphics("{fig1}")',
                    f'knitr::include_graphics("{fig2}")',
                    '```']
            for i in info:
                print(i, file=o_file, end="\n")

    print("\n\n", file=o_file)

def find_famfile(i_filelist, i_opt):
    if i_opt == 8:
        outpattern = []
        outfile = None

        for i in i_filelist:
            patt = re.search(r'\w+\.\w+_family_', i)
            if patt == None:
                pass
            else:
                outpattern.append(patt.group(0))
                outfile = i

        if len(outpattern) > 1:
            raise ValueError("Something is wrong, more than one file was found corresponding to 2 families!")
        elif len(outpattern) == 1:
            return outpattern[0], outfile
        else:
            return None, None
    elif i_opt == 7:
        outfile = None
        outpattern = set()
        
        for i in i_filelist:
            if "merged_family_GeneLocation" in i:
                outfile = i
            else:
                temp = os.path.basename(i)
                if "_family_" in temp:
                    FamID = temp.split("_family_")[0]
                    outpattern.add(FamID)
                
        if outfile != None:
            return list(outpattern), outfile
        else:
            raise ValueError("Unknown error")
    else:
        emsg = f"Unknown 'i_opt' value: {i_opt}"
        raise ValueError(emsg)
    
    
# Sort summary tables by name
def sort_tables(ilist, i_opt=None):
    # ilist = [i for i in ilist if "GeneOrganizationGenome" not in i]

    if i_opt == None:
        D = {0 : "", 1 : "", 2 : "", 3 : ""}
        Dtitles = {0 : "**Table 1. Genome-wide organization of gene family clusters**", \
                   1 : "**Table 2. Gene family cluster organization per scaffold**", \
                   2 : "**Table 3. Gene family cluster sizes**", \
                   3 : "**Table 4. Gene family membership**"}

        for yi in ilist:
            if "GeneOrganizationGenome" in yi:
                D[0] = yi

            elif "GeneOrganizationSummary" in yi:
                D[1] = yi

            elif "ClusterSizes" in yi:
                D[2] = yi

            elif "GeneLocation" in yi:
                D[3] = yi

        return D, Dtitles
    
    elif i_opt == "TwoFam":
        if len(ilist) == 8:
            TwoFam_hintname, TwoFamSumfile  = find_famfile(ilist, 8)

            if TwoFam_hintname == None: # checkpoint
                raise ValueError("Something is wrong, two family sum files not found")
            
            D = {0 : TwoFamSumfile, 1 : "", 2 : "",
                 3 : "", 4 : "", 5 : "", 6 : "", 7 : ""}
            Dtitles = {0 : "**Table 1. Two gene families cluster organization **", \
                       1 : "**Table 2a. Genome-wide organization of Gene family clusters**", \
                       2 : "**Table 2b. Gene family cluster organization per scaffold**", \
                       3 : "**Table 2c. Gene family cluster sizes **", \
                       4 : "**Table 3a. Genome-wide organization of Gene family clusters**", \
                       5 : "**Table 3b. Gene family cluster organization per scaffold**", \
                       6 : "**Table 3c. Gene family cluster sizes **", \
                       7 : "**Table 4. All gene families membership **"}
            
            # Get family names from hint
            temp = TwoFam_hintname.split("_family_")[0].split(".")
            fam1, fam2 = [*temp]
            
            for yi in ilist:
                if "GeneOrganizationG" in yi and os.path.basename(yi).startswith(f"{fam1}_"):
                    didx = 1 # dict index
                    D[didx] = yi
                    Dtitles[didx] = Dtitles[didx].replace("Gene family", f"{fam1} gene family")
                    
                elif "GeneOrganizationG" in yi and os.path.basename(yi).startswith(f"{fam2}_"):
                    didx = 4 # dict index
                    D[didx] = yi
                    Dtitles[didx] = Dtitles[didx].replace("Gene family", f"{fam2} gene family")

                elif "GeneOrganizationS" in yi and os.path.basename(yi).startswith(f"{fam1}_"):
                    didx = 2 # dict index
                    D[didx] = yi
                    Dtitles[didx] = Dtitles[didx].replace("Gene family", f"{fam1} gene family")
                    
                elif "GeneOrganizationS" in yi and os.path.basename(yi).startswith(f"{fam2}_"):
                    didx = 5 # dict index
                    D[didx] = yi
                    Dtitles[didx] = Dtitles[didx].replace("Gene family", f"{fam2} gene family")
                    
                elif "ClusterSizes" in yi and os.path.basename(yi).startswith(f"{fam1}_"):
                    didx = 3 # dict index
                    D[didx] = yi
                    Dtitles[didx] = Dtitles[didx].replace("Gene family", f"{fam1} gene family")

                elif "ClusterSizes" in yi and os.path.basename(yi).startswith(f"{fam2}_"):
                    didx = 6 # dict index
                    D[didx] = yi
                    Dtitles[didx] = Dtitles[didx].replace("Gene family", f"{fam2} gene family")
                    
                elif "GeneLocation" in yi:
                    didx = 7 # dict index
                    D[didx] = yi
                    Dtitles[didx] = Dtitles[didx].replace("All gene families", f"{fam1} and {fam2} gene families")

            return D, Dtitles
        
        elif len(ilist) == 7:
            Fams_hintname, FamSumfile  = find_famfile(ilist, 7)
            fam1, fam2 = [*Fams_hintname]
            
            D = {0 : "", 1 : "", 2 : "",
                 3 : "", 4 : "", 5 : "", 6 : ""}

            Dtitles = {0 : "**Table 1a. Genome-wide organization of Gene family clusters**", \
                       1 : "**Table 1b. Gene family cluster organization per scaffold**", \
                       2 : "**Table 1c. Gene family cluster sizes **", \
                       3 : "**Table 2a. Genome-wide organization of Gene family clusters**", \
                       4 : "**Table 2b. Gene family cluster organization per scaffold**", \
                       5 : "**Table 2c. Gene family cluster sizes **", \
                       6 : "**Table 3. All gene families membership **"}

            for yi in ilist:
                if "GeneOrganizationG" in yi and os.path.basename(yi).startswith(f"{fam1}_"):
                    didx = 0 # dict index
                    D[didx] = yi
                    Dtitles[didx] = Dtitles[didx].replace("Gene family", f"{fam1} gene family")
                    
                elif "GeneOrganizationG" in yi and os.path.basename(yi).startswith(f"{fam2}_"):
                    didx = 3 # dict index
                    D[didx] = yi
                    Dtitles[didx] = Dtitles[didx].replace("Gene family", f"{fam2} gene family")

                elif "GeneOrganizationS" in yi and os.path.basename(yi).startswith(f"{fam1}_"):
                    didx = 1 # dict index
                    D[didx] = yi
                    Dtitles[didx] = Dtitles[didx].replace("Gene family", f"{fam1} gene family")
                    
                elif "GeneOrganizationS" in yi and os.path.basename(yi).startswith(f"{fam2}_"):
                    didx = 4 # dict index
                    D[didx] = yi
                    Dtitles[didx] = Dtitles[didx].replace("Gene family", f"{fam2} gene family")
                    
                elif "ClusterSizes" in yi and os.path.basename(yi).startswith(f"{fam1}_"):
                    didx = 2 # dict index
                    D[didx] = yi
                    Dtitles[didx] = Dtitles[didx].replace("Gene family", f"{fam1} gene family")

                elif "ClusterSizes" in yi and os.path.basename(yi).startswith(f"{fam2}_"):
                    didx = 5 # dict index
                    D[didx] = yi
                    Dtitles[didx] = Dtitles[didx].replace("Gene family", f"{fam2} gene family")
                    
                elif "GeneLocation" in yi:
                    didx = 6 # dict index
                    D[didx] = yi
                    Dtitles[didx] = Dtitles[didx].replace("All gene families", f"{fam1} and {fam2} gene families")
                    
            return D, Dtitles
        else:
            raise ValueError("Something is wrong, unknown error")
    
    else:
        raise ValueError("Unknown 'i_opt'")

def Sort_IndividualPlots(i_list, i_guide, i_hint, i_cutname=False):
    tempDF = pd.DataFrame()

    if i_list[0].startswith(i_hint):
        shortname = []
        for i in i_list:
            shortname.append(i.replace("merged_","").replace(".svg",""))

        tempDF["rawnames"] = i_list
        tempDF["shortnames"] = shortname

        # Create the dictionary that defines the order for sorting
        sorterIndex = dict(zip(i_guide, range(len(i_guide))))

        # Generate a rank column that will be used to sort
        # the dataframe numerically
        tempDF['order_indx'] = tempDF['shortnames'].map(sorterIndex)
        tempDF = tempDF.sort_values(by="order_indx")
    
    elif i_hint == "_2fam_" or i_hint == "_fam_":
        if i_hint in i_list[0]:
            shortname = []
            for i in i_list:
                shortname.append(i.split(i_hint)[-1].replace(".svg",""))
            # print(shortname)
            tempDF["rawnames"] = i_list
            tempDF["shortnames"] = shortname

            # Create the dictionary that defines the order for sorting
            sorterIndex = dict(zip(i_guide, range(len(i_guide))))

            # Generate a rank column that will be used to sort
            # the dataframe numerically
            tempDF['order_indx'] = tempDF['shortnames'].map(sorterIndex)
            tempDF = tempDF.sort_values(by="order_indx")
    
    elif i_hint == "Singlefamily":
        shortname = []
        for i in i_list:
            shortname.append(os.path.basename(i).replace(".svg",""))
        
        if i_cutname == True:
            shortname = [os.path.basename(i).split("_matrices.")[-1].split(".matrix.")[0] for i in shortname]
        tempDF["rawnames"] = i_list
        tempDF["shortnames"] = shortname

        # Create the dictionary that defines the order for sorting
        sorterIndex = dict(zip(i_guide, range(len(i_guide))))

        # Generate a rank column that will be used to sort
        # the dataframe numerically
        tempDF['order_indx'] = tempDF['shortnames'].map(sorterIndex)
        tempDF = tempDF.sort_values(by="order_indx")        
    else:
        raise ValueError("Unknown 'i_hint' value")
        
    sortednames = tempDF["rawnames"].to_list()
    return sortednames


# In[ ]:


# Fx: for two families analysis

def sort_2FamClustSum_plots(i_obj, i_opt):

    if i_opt == "SummaryPlots":
        sorted_list = []
        famnames = []
        CasesTitles = []
        
        # Find the hints for "merged analysis bw 2 families"
        for i in i_obj:
            iname = os.path.basename(i)
            if iname == "merged":
                sorted_list.append(i)
                CasesTitles.append("A. Two families merged analysis")
                
        # Find the hints for "interaction analysis bw 2 families"
        for i in i_obj:
            iname = os.path.basename(i)
            if "." in iname:
                sorted_list.append(i)
                famnames = iname.split(".")
                CasesTitles.append("B. Two families joint analysis")

            # if there's no interaction, only plots for merged analysis and individual families will be printed in the final report
            if len(famnames) == 0:
                nointeraction = True
            else:
                nointeraction = False
        
        # if there's no interaction, find the family IDs, using the basename of other plots
        if nointeraction == True:
            for i in i_obj:
                iname = os.path.basename(i)
                if iname == "merged":
                    pass
                elif "." in iname:
                    pass
                else:
                    famid = os.path.basename(i)
                    if famid not in famnames: famnames.append(famid)

        for i in i_obj:
            iname = os.path.basename(i)
            if iname == "merged":
                pass
            elif "." in iname:
                pass
            else:
                famnames = sorted(famnames)
                for idx, ifam in enumerate(famnames):
                    if iname.endswith(ifam):
                        sorted_list.append(i)
                        if idx == 0:
                            if nointeraction == False:
                                CasesTitles.append(f"C. {ifam} family clusters")
                            elif nointeraction == True:
                                CasesTitles.append(f"B. {ifam} family clusters")
                        elif idx == 1:
                            if nointeraction == False:
                                CasesTitles.append(f"C. {ifam} family clusters")
                            elif nointeraction == True:
                                CasesTitles.append(f"B. {ifam} family clusters")
                        else:
                            raise ValueError("Something is wrong, only two families are expected")
        return sorted_list, CasesTitles
    
    elif i_opt == "IndividualPlots":
        sorted_list = []
        famnames = []
        CasesTitles = []
        
        for i in i_obj:
            iname = os.path.basename(i)
            if iname == "merged":
                sorted_list.append(i)
                CasesTitles.append("A. Two families merged analysis")
                
        for i in i_obj:
            iname = os.path.basename(i)
            if "." in iname:
                sorted_list.append(i)
                famnames = iname.split(".")
                CasesTitles.append("B. Two families joint analysis")
                
            # if there's no interaction, only plots for merged analysis and individual families will be printed in the final report
            if len(famnames) == 0:
                nointeraction = True
            else:
                nointeraction = False
        
        # if there's no interaction, find the family IDs, using the basename of other plots
        if nointeraction == True:
            for i in i_obj:
                iname = os.path.basename(i)
                if iname == "merged":
                    pass
                elif "." in iname:
                    pass
                else:
                    famid = os.path.basename(i)
                    if famid not in famnames: famnames.append(famid)

        for i in i_obj:

            iname = os.path.basename(i)
            if iname == "merged":
                pass
            elif "." in iname:
                pass
            else:
                famnames = sorted(famnames)
                for idx, ifam in enumerate(famnames):
                    if iname.endswith(ifam):
                        sorted_list.append(i)
                        if idx == 0:
                            if nointeraction == False:
                                CasesTitles.append(f"C. {ifam} family clusters")
                            elif nointeraction == True:
                                CasesTitles.append(f"B. {ifam} family clusters")
                        elif idx == 1:
                            if nointeraction == False:
                                CasesTitles.append(f"C. {ifam} family clusters")
                            elif nointeraction == True:
                                CasesTitles.append(f"B. {ifam} family clusters")
                        else:
                            raise ValueError("Something is wrong, only two families are expected")                        
        return sorted_list, CasesTitles
    
    else:
        raise ValueError("Unknown 'i_opt' value")
    

def custom_sorter(i_element):
    if "Global_ClusterSize" in i_element:
        return 0
    elif "ClusterSize_Distribution" in i_element:
        return 1
    else:
        return 2

def write_figures_merged_as_singlefigs(i_dict, i_key, o_file, i_echo, i_opt):
    if i_echo == False:
        i_echo = "FALSE"
    elif i_echo == True:
        i_echo = "TRUE"
        
    temp, titles = sort_2FamClustSum_plots(i_dict[i_key], i_opt)
    D_titles = {k : v for k,v in zip(temp, titles)}
    
    for idir in temp:
        print(f"#### **{D_titles[idir]}**", file=o_file, end="\n")
        reordered_filelist = list(os.listdir(idir))[::-1]
        if i_opt == "IndividualPlots":
            if reordered_filelist[0].startswith("merged_"):
                hint = "merged_"
                reordered_filelist = Sort_IndividualPlots(reordered_filelist, guidelist, hint)
            elif "_2fam_" in reordered_filelist[0]:
                hint = "_2fam_"
                reordered_filelist = Sort_IndividualPlots(reordered_filelist, guidelist, hint)
            elif "_fam_" in reordered_filelist[0]:
                hint = "_fam_"
                reordered_filelist = Sort_IndividualPlots(reordered_filelist, guidelist, hint)
            else:
                raise ValueError("Unknown error")
            
            for idx, fig in enumerate(reordered_filelist):
                fig = f"{idir}/{fig}"

                info = [f'**Figure {idx+1}: {os.path.basename(fig)}**', '<div>\n',
                            '```{' + 'r, out.width="75%", echo={}'.format(i_echo) +'}',
                            'library(knitr)',
                            f'knitr::include_graphics("{fig}")',
                            '```\n\n']
                for i in info:
                    print(i, file=o_file, end="\n")

        elif i_opt == "SummaryPlots":
            reordered_filelist = sorted(reordered_filelist, key=custom_sorter) # Sort the output graphs by their name, first GlobalClusterSize, then ClusteSize_distribution

            for idx, fig in enumerate(reordered_filelist):

                fig = f"{idir}/{fig}"
                info = [f'**Figure {idx+1}: {os.path.basename(fig)}**', '<div>\n',
                        '```{' + 'r, out.width="75%", echo={}'.format(i_echo) +'}',
                        'library(knitr)',
                        f'knitr::include_graphics("{fig}")',
                        '```\n\n']

                for i in info:
                    print(i, file=o_file, end="\n")
        else:
            raise ValueError("Unknown error")
                
# Fx to create a report 
def CreateReport(i_dict, iFamily, igval, o_file, i_opt=None):
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Header
    ReportFile = open(o_file, "w")
    
    maintitle = f"{iFamily}_{igval}"

    header_info = ["---", "title: GALEON", 
                   "output: html_document", 
                   "---", f"## Report: {maintitle} " + "{.tabset}"]
    for i in header_info:
        print(i, file=ReportFile, end="\n")
            
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Open a new worksheet
    temp = "Summary Tables"
    WorkSheetName = f"### {temp}"

    print(WorkSheetName, file=ReportFile, end="\n")
    
    # Create a table
    d_tables, d_tablenames = sort_tables(i_dict["SummaryFiles"], i_opt)

    for (ikey, table), (jkey, table_name) in zip(d_tables.items(), d_tablenames.items()):
        print(table_name, file=ReportFile, end="\n")
        
        if ikey == 0:
            write_table(table, ReportFile, 1, echo_paths)
        else:
            write_table(table, ReportFile, 2, echo_paths)


    if MannWhitney_results == True:
        print("**Table 5. Cst index at genome level**", file=ReportFile, end="\n")
        write_table(i_dict["Cst glob stats"], ReportFile, 2, echo_paths)

        print("**Table 6. Mann-Whitney test results**", file=ReportFile, end="\n")
        write_table(i_dict["MW results"], ReportFile, 2, echo_paths)
        
        print("**Table 7. Mann-Whitney test raw data**", file=ReportFile, end="\n")
        write_table(i_dict["MW raw data"], ReportFile, 2, echo_paths)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Open a new worksheet
    temp = "Summary Plots"
    WorkSheetName = f"### {temp}"

    print(WorkSheetName, file=ReportFile, end="\n")

    # Export figures
    # write_figures_merged_as_twofigs(i_dict, f"SummaryPlots_{igval}", ReportFile)
    if i_opt == "TwoFam":
        write_figures_merged_as_singlefigs(i_dict, f"SummaryPlots_{igval}", ReportFile, echo_paths, "SummaryPlots")
    else:
        write_figures(i_dict, f"SummaryPlots_{igval}", ReportFile, echo_paths, "SummaryPlots")
        
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Open a new worksheet
    temp = "Individual Plots"
    WorkSheetName = f"### {temp}"

    print(WorkSheetName, file=ReportFile, end="\n")

    # Export figures
    if i_opt == "TwoFam":
        write_figures_merged_as_singlefigs(i_dict, f"IndividualPlots_{igval}", ReportFile, echo_paths, "IndividualPlots")
    else:
        write_figures(i_dict, f"IndividualPlots_{igval}", ReportFile, echo_paths, "IndividualPlots")
        
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Open a new worksheet
    temp = "Physical distance heatmaps"
    WorkSheetName = f"### {temp}"

    print(WorkSheetName, file=ReportFile, end="\n")

    # Export figures
    if i_opt == "TwoFam":
        write_PhysMxHeatmap_figures(i_dict, "GaleonClusterFinderResults", ReportFile, echo_paths)
    else:
        write_PhysMxHeatmap_figures(i_dict, "GaleonClusterFinderResults", ReportFile, echo_paths, False)
    
    
    if EvoDistUse == True:
        # Open a new worksheet
        temp = "Physical and Evolutionary plots"
        WorkSheetName = f"### {temp}"

        print(WorkSheetName, file=ReportFile, end="\n")

        # Export figures
        temp2 = list(i_dict["PhysicalEvoDistance_Figures"].keys())
        scf_order = Sort_IndividualPlots(temp2, scfguidelist, "Singlefamily", True)
        
        reord_dict = {key : i_dict["PhysicalEvoDistance_Figures"][key] for key in scf_order}
        
        write_2_figures(reord_dict, iFamily, igval, echo_paths, ReportFile)   
        
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Open a new worksheet
    temp = "HELP"
    WorkSheetName = f"### {temp}"

    print(WorkSheetName, file=ReportFile, end="\n")

    with open(f"{PATHtoGaleonScripts}/Help.Rmd") as helpfile:
        for line in helpfile:
            if "Help/Figures" in line:
                line = line.replace("Help/Figures", f"{PATHtoGaleonScripts}/Help/Figures")
                ReportFile.write(line)
            else:
                ReportFile.write(line)

    ReportFile.close()
    


# In[ ]:


# Create reports in Rmn format for each family
for Family, v in D_GLOB.items():
    for gval, info in v.items():
        reportname = f"{Family}_{gval}_Report.Rmd"
        print(Family, gval)
        print(f"- {reportname}")
        
        if Family.startswith("merged_"):
            CreateReport(info, Family, gval, reportname, "TwoFam")
        else:
            CreateReport(info, Family, gval, reportname)
        
        renderscript_name = f"{Family}_{gval}_render_script.R"
        with open(renderscript_name, "w") as o1:
            o1.write(f'rmarkdown::render("{reportname}", "html_document")')


# Render in html format
for script in [i for i in os.listdir(".") if i.endswith("render_script.R")]:
    subprocess.call(["Rscript", script])
    os.remove(script)



# Move to the output directory
for htmlfile in [i for i in os.listdir(".") if i.endswith(".html")]:
    shutil.move(htmlfile, OUTdir)

# Move to the output directory
for Rmdfile in [i for i in os.listdir(".") if i.endswith("Report.Rmd")]:
    shutil.move(Rmdfile, OUTdir)

# Move the report dir to the main results dir
if os.path.exists(f"{GaleonClusterFinderResults}/{OUTdir}"):
    shutil.rmtree(f"{GaleonClusterFinderResults}/{OUTdir}")
    shutil.move(OUTdir, GaleonClusterFinderResults)
else:
    shutil.move(OUTdir, GaleonClusterFinderResults)
