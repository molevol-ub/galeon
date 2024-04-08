import pandas as pd
import seaborn as sns
import numpy as np
import os, sys
import matplotlib.pyplot as plt
from collections import Counter


# Input data 
GeneFamilyFile = sys.argv[1] # "OutTable.GR.tsv"
hint = sys.argv[2]
PlotDir = sys.argv[3]
# (Optional) Retain only the information of the main scaffolds (chromosomes)
ScfToRetain = sys.argv[4]
ColorOpt = sys.argv[5]
ColorValueChoice = sys.argv[6]

FAMname = os.path.basename(GeneFamilyFile).split("_")[0]
FAMname

# Column containing the IDs of scaffolds
ScaffoldIDcol = "ScfName"


# checkpoint
if ScaffoldIDcol not in ["ScfName", "ScaffoldID"]:
    emsg = f"Unknown name for 'ScaffoldIDcol={ScaffoldIDcol} column. Accepted values: ['ScfName', 'ScaffoldID'] "
    raise ValueError(emsg)


# Output dirs
SummaryPlotsDir = f"{PlotDir}/SummaryPlots_{hint}"
SinglePlotsDir = f"{PlotDir}/IndividualPlots_{hint}"

for i in [SummaryPlotsDir, SinglePlotsDir]:
    if os.path.exists(i):
        pass
    else:
        os.mkdir(i)


# In[ ]:


# Check Scaffold to retain 
if ScfToRetain == "ALL":
    print(f"'ScfToRetain' option 'ALL': {ScfToRetain} - all scaffolds are considered")
elif os.path.exists(ScfToRetain):
    print(f"'ScfToRetain' option 'FILE': {ScfToRetain} - all scaffolds from the input file are considered")
elif isinstance(ScfToRetain, str) == True:
    try:
        ScfToRetain = int(ScfToRetain)
        if isinstance(ScfToRetain, int):
            print(f"'ScfToRetain' option 'NUM': {ScfToRetain} - only the first {ScfToRetain} largest scaffolds are considered")
    except ValueError as e:
        emsg = f"Error in 'ScfToRetain', unknown input type: {type(ScfToRetain)}. Allowed input: \n\t\
            -'ALL' - all scaffolds are considered \n\t\
            -'FILE' - all scaffolds from the input file are considered\n\t\
            -'NUM' - only the first NUM largest scaffolds are considered' \n"
        raise ValueError(emsg)
    


# In[ ]:


# Load dataframe
DF = pd.read_csv(GeneFamilyFile, sep="\t")

# Fx to filter retain only some scaffolds of interest 
def GetSomeScf(i_DF, iScfIDcolumn, i_filter):
    if isinstance(i_filter, int):

        D_scflen = {} # Save scf lengths in this dictionary

        for k, v in i_DF.groupby([iScfIDcolumn,"ScfLength"]):
            D_scflen[k[0]] = k[1]

        # Sort this dict by scaffold length
        D_scflen_sorted = sorted(D_scflen.items(), key=lambda x:x[1], reverse=True)

        # Extract the 'x' scaffolds of interest
        ScfOfInterest = [i[0] for i in D_scflen_sorted[:i_filter]]

        # Export only a part of dataframe that contains the scf of interest
        oDF = i_DF[i_DF[iScfIDcolumn].isin(ScfOfInterest)].copy(deep=True)
        return oDF
    
    else:
        with open(ScfToRetain) as f1:
            ScfList2retain = f1.read().strip().split("\n")
            
            oDF = i_DF[i_DF["ScaffoldID"].isin(ScfList2retain)]
            
            return oDF

if ScfToRetain == "ALL":
    DFcp = DF.copy(deep=True)
else:
    DFcp = GetSomeScf(DF, ScaffoldIDcol, ScfToRetain)


# In[ ]:


def GroupBy_ChrClusterID(i_DF, iScfIDcolumn):
    # Group by Chromosome and ClusterID
    clustersizes = []

    # Create records by chr
    outrecords = []

    DFby = i_DF.groupby([iScfIDcolumn, "ClusterID"])
    for k,v in DFby:
        genes_num = len(v["GeneID"])
        clustersizes.append(genes_num)
        # print(k, genes_num)
        temp = list(k) + [genes_num]
        outrecords.append(temp)

    # Create a dataframe with chr - clusterID - gene num info

    DF_Chr = pd.DataFrame.from_records(outrecords)

    DF_Chr.columns = [iScfIDcolumn, "ClusteredGenesNumber", "GeneNumber"]
    return DF_Chr, clustersizes
    
# Create a dataframe with the following information:
# Number of genes per cluster in each chr/scf
DF_Chr, ClusterSizes = GroupBy_ChrClusterID(DFcp, ScaffoldIDcol)


# In[5]:


def GeneralSummary_Plot(iClusterSizes, iColor, iPlotName):
    
    # Counter clusters size frequency: How many cluster there are of 'x' size? Ex: 20 clusters of 5 genes, 10 clusters of 4 genes
    D_Count = dict(Counter(iClusterSizes))

    # Create a dataframe
    outrecords = []
    for k,v in D_Count.items():
        outrecords.append([k,v])
    
    # Create a DF with frequency of cluster sizes
    DF_byFreqGenenum = pd.DataFrame.from_records(outrecords)
    DF_byFreqGenenum.columns = [f"Length of {FAMname} array", "Frequency"]
    DF_byFreqGenenum

    # Create a Barplot
    fig, ax1 = plt.subplots(figsize=(12, 8))
    BarPlot = sns.barplot(data=DF_byFreqGenenum, 
                          x = f"Length of {FAMname} array", 
                          y = "Frequency", 
                          edgecolor="black",
                          color=iColor,
                          ax=ax1, 
                          errorbar=None,
                          width=0.8)
    
    # Set y ticks
    new_y_ticks = range(1, max(DF_byFreqGenenum["Frequency"])+1, 5)  # Replace this list with your desired tick positions
    plt.yticks(new_y_ticks)
    plt.xlabel(f'Length of {FAMname} array', fontsize=16)
    plt.ylabel('Frequency', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)

    # Add freq as lables on the top of the BarPlot
    for p in BarPlot.patches:
        BarPlot.annotate(int(p.get_height()), (p.get_x() + p.get_width() / 2., p.get_height()),  
                           ha = 'center', 
                           va = 'center', 
                           xytext = (0, 10), 
                           textcoords = 'offset points',
                           fontsize = 12)
    ax1.spines[['right', 'top']].set_visible(False)
    plt.savefig(f"{SummaryPlotsDir}/{iPlotName}")
    # plt.show()
    # plt.close()
    plt.clf()

# Create a General Bar Plot (Frequency of cluster sizes)
# GeneralSummary_Plot(ClusterSizes, "lightgrey", "Global_ClusterSize_Distribution.png")
GeneralSummary_Plot(ClusterSizes, "lightgrey", f"{FAMname}_Global_ClusterSize_Distribution.svg")


# In[6]:


def ModDf(i_DF, iScfIDcolumn, iNumOFClusters):
    # Add scaffolds as keys
    D_temp = {k : {} for k in i_DF[iScfIDcolumn].unique()} # temp dictionary. It will be used to create a barplot with cluster frequencies for each chromosome

    # Add Cluster sizes as subkeys
    for k,v in D_temp.items():
        for kk in iNumOFClusters:
            v[kk] = 0

    for idx, r in i_DF.iterrows():
        scfID, clustSize = r[iScfIDcolumn], r["ClusterSize"]
        D_temp[scfID][clustSize] = r["Frequency"]

    # Dict to dataframe
    # Convert the dictionary to a list of dictionaries
    data = [{iScfIDcolumn: chromosome, "ClusterSize": position, "Frequency": value}
            for chromosome, positions in D_temp.items()
            for position, value in positions.items()]

    # Create a DataFrame
    df = pd.DataFrame(data)

    return df



# In[7]:


def GeneralBar_Plot(iDF_Chr, iScaffoldIDcolumn, oDirName, i_color_opt, i_colorOBJ, iVerticalLine, iZeroLineStyle, iZeroLineSize, iPlotName, izerofreq_plot_opt=True):
    # Create a dataframe grouping by Chromosome and GeneNumber
    DF_byChrGeneNum = iDF_Chr.groupby([iScaffoldIDcolumn, "GeneNumber"])

    # Create another DF with the frequency of gene cluster sizes per chromosome

    outrecords = []
    for k,v in DF_byChrGeneNum:
        # print(k)
        temp = list(k) + [v.shape[0]]
        outrecords.append(temp)

    DF_byChrFreqGenenum = pd.DataFrame.from_records(outrecords)
    DF_byChrFreqGenenum.columns = [iScaffoldIDcolumn, "ClusterSize", "Frequency"]
    
    # Num of clusters categories (each category is a cluster size)
    NumOFClusters = sorted(set(DF_byChrFreqGenenum["ClusterSize"].to_list()))
    
    # Create a new dataframe, with 0 frequencies (with plotting purpose only)
    newDF = ModDf(DF_byChrFreqGenenum, iScaffoldIDcolumn, NumOFClusters)
    
    # Create a Barplot
    MultipleScf_Barplot(DF_byChrFreqGenenum, iScaffoldIDcolumn, newDF, NumOFClusters, oDirName, i_color_opt, i_colorOBJ, iVerticalLine, iZeroLineStyle, iZeroLineSize, iPlotName, izerofreq_plot_opt)
    
    return DF_byChrFreqGenenum, NumOFClusters

def MultipleScf_Barplot(iFreqDF, iScafColName, inewDF, inumclust, oDir, icolor_opt, icolorOBJ, iVertical, iZeroLineStyle, iZeroLineSize, iplotName, zero_plot_opt=True):
    
    # Create a barplot with several scaffolds plotted on the same "cluster size" category
    fig, ax1 = plt.subplots(figsize=(12, 8))
    
    inewDF.columns = ["ScfName", f"Length of {FAMname} array", "Frequency"]
    if icolor_opt == "OneColor":
        BarPlot = sns.barplot(data=inewDF, 
                              x = f"Length of {FAMname} array", 
                              y = "Frequency", 
                              hue = iScafColName,
                              color=icolorOBJ,
                              edgecolor="black",
                              ax=ax1, 
                              errorbar=None,
                              width=0.8)


    elif icolor_opt == "Palette":
        BarPlot = sns.barplot(data=inewDF, 
                              x = f"Length of {FAMname} array", 
                              y = "Frequency", 
                              hue = iScafColName,
                              palette=icolorOBJ,
                              edgecolor="black",
                              ax=ax1, 
                              errorbar=None,
                              width=0.8)

    else:
        emsg = f"Unknown 'i_color_opt'={icolor_opt} in 'MultipleScf_Barplot' function. Admitted values: ['OneColor', 'Palette']"
        raise ValueError(emsg)
        
    # Set ticks on y axis
    new_y_ticks = range(1, max(iFreqDF["Frequency"])+1)  # Replace this list with your desired tick positions
    plt.yticks(new_y_ticks)
    plt.xlabel(f'Length of {FAMname} array', fontsize=16)
    plt.ylabel('Frequency', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.legend(title='Scaffold ID', fontsize=14)

    # Add labels (Frequency values) on the top of bars
    for p in BarPlot.patches:
        if zero_plot_opt == True:
            if np.isnan(p.get_height()):
                pass
            else:
                BarPlot.annotate(int(p.get_height()), (p.get_x() + p.get_width() / 2., p.get_height()),  
                                   ha = 'center', 
                                   va = 'center', 
                                   xytext = (0, 10), 
                                   textcoords = 'offset points',
                                   fontsize = 12)
        elif zero_plot_opt == False:
            if np.isnan(p.get_height()):
                pass
            elif p.get_height() == 0:
                pass
            else:
                BarPlot.annotate(int(p.get_height()), (p.get_x() + p.get_width() / 2., p.get_height()),  
                                   ha = 'center', 
                                   va = 'center', 
                                   xytext = (0, 10), 
                                   textcoords = 'offset points',
                                   fontsize = 12)
    
    
    for ix, a in enumerate(ax1.patches):
        if inewDF.loc[ix, 'Frequency'] == 0:
            x_start = a.get_x()
            width = a.get_width()

            ax1.plot([x_start, x_start+width], 2*[inewDF.loc[ix, 'Frequency']], 
                     iZeroLineStyle, 
                     linewidth=iZeroLineSize, c='k')
    
    
    # Add vertical separators on x axis
    if iVertical == "AddVerticalLine":
        for i in range(len(inumclust)-1):
            ax1.axvline(i+0.5, color='gray', linestyle='-', linewidth=0.5)
        ax1.spines[['right', 'top']].set_visible(False)
        plt.savefig(f"{oDir}/{iplotName}")
        # plt.show()
        # plt.close()
        plt.clf()

    elif iVertical == "False":
        ax1.spines[['right', 'top']].set_visible(False)
        plt.savefig(f"{oDir}/{iplotName}")        
        # plt.show()
        # plt.close()
        plt.clf()
    
    else:
        emsg = f"Unknown 'iVertical'={iVertical} in 'MultipleScf_Barplot' function. Admitted values: ['AddVerticalLine', 'False']"
        raise ValueError(emsg)        


    
    
# DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, SummaryPlotsDir, "OneColor", "red", "False", "--", 5, "ClusterSize_Distribution_by_Scaffolds.png", False)
# Change of visualization of zero frequency lines depending on the number of chromosomes (with visualization purpose only)
if len(DF_Chr["ScfName"].unique()) < 3:
    DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, SummaryPlotsDir, ColorOpt, ColorValueChoice, "False", "-", 5, f"{FAMname}_ClusterSize_Distribution_by_Scaffolds.svg", False)
else:
    DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, SummaryPlotsDir, ColorOpt, ColorValueChoice, "False", "--", 5, f"{FAMname}_ClusterSize_Distribution_by_Scaffolds.svg", False)
# DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, SummaryPlotsDir, "Palette", "husl", "AddVerticalLine", "--", 5, "test2.png", False)


# In[ ]:





# In[8]:


def ScfbyScf_Barplot(iDF_byChrFreqGenenum, iScaffoldIDcolname, iNumOfClusters, oDirname, icolor_opt, icolorOBJ, iVertical, iZeroLineStyle, iZeroLineSize):
    DF_byChrFreqGenenum_gb = iDF_byChrFreqGenenum.groupby(iScaffoldIDcolname)

    for kScf, vDF in DF_byChrFreqGenenum_gb:
        # Out plot name and format
        Plotname = kScf + ".svg"
        
        # Create a new dataframe, with 0 frequencies (with plotting purpose only)
        newDF = ModDf(vDF, iScaffoldIDcolname, iNumOfClusters)

        # Create a Barplot
        MultipleScf_Barplot(vDF, iScaffoldIDcolname, newDF, iNumOfClusters, oDirname, icolor_opt, icolorOBJ, iVertical, iZeroLineStyle, iZeroLineSize, Plotname)
        
        
ScfbyScf_Barplot(DF_ChrFreqGenenum, ScaffoldIDcol, NumOfClusters, SinglePlotsDir, "OneColor", "lightgrey", "AddVerticalLine", "-", 5)

