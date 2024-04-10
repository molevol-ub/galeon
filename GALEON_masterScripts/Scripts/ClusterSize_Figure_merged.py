import pandas as pd
import seaborn as sns
import numpy as np
import os, sys
import matplotlib.pyplot as plt
from collections import Counter

GeneFamilyFile = sys.argv[1] # "OutTable.GR.tsv"
hint = sys.argv[2]
PlotDir = sys.argv[3]

# (Optional) Retain only the information of the main scaffolds (chromosomes)
ScfToRetain = sys.argv[4]
ColorOpt = sys.argv[5]
ColorValueChoice = sys.argv[6]

# GeneFamilyFile = "Plots/merged_fam/merged_family_GeneLocation.1.0g.tsv" # "OutTable.GR.tsv"
# hint = "1.0g"
# PlotDir = "Plots/merged_fam"
# ScfToRetain = "7"

FAMname = os.path.basename(GeneFamilyFile).split("_")[0]

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
    

# Load dataframe
DF = pd.read_csv(GeneFamilyFile, sep="\t")


def SplitDF_by_Fam(iDF):
    # Get the list of clusters composed of Two families
    TwoFams_ClustersEntries = []
    for k, vDF in iDF.groupby(["ScaffoldID", "ClusterID"]):
        if len(vDF["FamID"].unique()) != 1:
            # print(k, vDF["FamID"].unique())
            TwoFams_ClustersEntries.append(list(k))


    # Split the input DataFrame into one corresponding to clusters composed of two families and 
    # one corresponding to clusters composed of single families

    outlist_2fam = []
    outlist_1fam = []

    for idx, row in iDF.iterrows():
        if [row["ScaffoldID"], row["ClusterID"]] in TwoFams_ClustersEntries:
            outlist_2fam.append(row)
        else:
            outlist_1fam.append(row)

    oDF_1Fam = pd.DataFrame.from_records(outlist_1fam)
    oDF_2Fam = pd.DataFrame.from_records(outlist_2fam)
    
    # checkpoint
    if iDF.shape[0]- oDF_2Fam.shape[0] - oDF_1Fam.shape[0] != 0:
        raise ValueError("Something is wrong. The sum of clustered+unclustered dataframe size doesn't match")
    
    return oDF_1Fam, oDF_2Fam

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
    DF_1Fam, DF_2Fam = SplitDF_by_Fam(DF)
else:
    DFcp = GetSomeScf(DF, ScaffoldIDcol, ScfToRetain)
    DF_1Fam, DF_2Fam = SplitDF_by_Fam(DFcp)

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
### DF_Chr, ClusterSizes = GroupBy_ChrClusterID(DFcp, ScaffoldIDcol)


def GeneralSummary_Plot(iClusterSizes, iColor, iPlotName, iFAMname):
    
    # Counter clusters size frequency: How many cluster there are of 'x' size? Ex: 20 clusters of 5 genes, 10 clusters of 4 genes
    D_Count = dict(Counter(iClusterSizes))

    # Create a dataframe
    outrecords = []
    for k,v in D_Count.items():
        outrecords.append([k,v])
    
    # Create a DF with frequency of cluster sizes
    DF_byFreqGenenum = pd.DataFrame.from_records(outrecords)
    DF_byFreqGenenum.columns = [f"Length of {iFAMname} array", "Frequency"]
    DF_byFreqGenenum

    # Create a Barplot
    fig, ax1 = plt.subplots(figsize=(12, 8))
    BarPlot = sns.barplot(data=DF_byFreqGenenum, 
                          x = f"Length of {iFAMname} array", 
                          y = "Frequency", 
                          edgecolor="black",
                          color=iColor,
                          ax=ax1, 
                          errorbar=None,
                          width=0.8)
        
    if iFAMname == "merged":
        plt.xlabel('Length of two merged families array', fontsize=16)
        plt.ylabel('Frequency', fontsize=16)
        plt.tick_params(axis='both', which='major', labelsize=14)
    else:
        plt.xlabel(f'Length of {iFAMname} array', fontsize=16)
        plt.ylabel('Frequency', fontsize=16)
        plt.tick_params(axis='both', which='major', labelsize=14)

    # Set y ticks
    new_y_ticks = range(1, max(DF_byFreqGenenum["Frequency"])+1, 5)  # Replace this list with your desired tick positions
    plt.yticks(new_y_ticks)
    
    # Add freq as lables on the top of the BarPlot
    for p in BarPlot.patches:
        BarPlot.annotate(int(p.get_height()), (p.get_x() + p.get_width() / 2., p.get_height()),  
                           ha = 'center', 
                           va = 'center', 
                           xytext = (0, 10), 
                           textcoords = 'offset points',
                           fontsize = 12)
    
    # remove plot lines: right and top
    ax1.spines[['right', 'top']].set_visible(False)

    plt.savefig(f"{SummaryPlotsDir}/{iPlotName}")
    # plt.show()
    # plt.close()
    plt.clf()

    
# Create a General Bar Plot (Frequency of cluster sizes)
# GeneralSummary_Plot(ClusterSizes, "lightgrey", "Global_ClusterSize_Distribution.png")
### GeneralSummary_Plot(ClusterSizes, "lightgrey", "Global_ClusterSize_Distribution.svg")


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


def GeneralBar_Plot(iDF_Chr, iScaffoldIDcolumn, oDirName, i_color_opt, i_colorOBJ, iVerticalLine, iZeroLineStyle, iZeroLineSize, iPlotName, iFamName, izerofreq_plot_opt=True):
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
    MultipleScf_Barplot(DF_byChrFreqGenenum, iScaffoldIDcolumn, newDF, NumOFClusters, oDirName, i_color_opt, i_colorOBJ, iVerticalLine, iZeroLineStyle, iZeroLineSize, iFamName, iPlotName, izerofreq_plot_opt)
    
    return DF_byChrFreqGenenum, NumOFClusters


def MultipleScf_Barplot(iFreqDF, iScafColName, inewDF, inumclust, oDir, icolor_opt, icolorOBJ, iVertical, iZeroLineStyle, iZeroLineSize, iFamName, iplotName, zero_plot_opt=True):
    
    # Create a barplot with several scaffolds plotted on the same "cluster size" category
    fig, ax1 = plt.subplots(figsize=(12, 8))
    
    if icolor_opt == "OneColor":
        BarPlot = sns.barplot(data=inewDF, 
                              x = "ClusterSize", 
                              y = "Frequency", 
                              hue = iScafColName,
                              color=icolorOBJ,
                              edgecolor="black",
                              ax=ax1, 
                              errorbar=None,
                              width=0.8)
        if iFamName == "merged":
            plt.xlabel('Length of two merged families array', fontsize=16)
            plt.ylabel('Frequency', fontsize=16)
            plt.tick_params(axis='both', which='major', labelsize=14)
            plt.legend(title='Scaffold ID', fontsize=14)
        else:
            plt.xlabel(f'Length of {iFamName} array', fontsize=16)
            plt.ylabel('Frequency', fontsize=16)
            plt.tick_params(axis='both', which='major', labelsize=14)
            plt.legend(title='Scaffold ID', fontsize=14)

        

    elif icolor_opt == "Palette":
        BarPlot = sns.barplot(data=inewDF, 
                              x = "ClusterSize", 
                              y = "Frequency", 
                              hue = iScafColName,
                              palette=icolorOBJ,
                              edgecolor="black",
                              ax=ax1, 
                              errorbar=None,
                              width=0.8)
        
        if iFamName == "merged":
            plt.xlabel('Length of two merged families array', fontsize=16)
            plt.ylabel('Frequency', fontsize=16)
            plt.tick_params(axis='both', which='major', labelsize=14)
            plt.legend(title='Scaffold ID', fontsize=14)


        else:
            plt.xlabel(f'Length of {iFamName} array', fontsize=16)
            plt.ylabel('Frequency', fontsize=16)
            plt.tick_params(axis='both', which='major', labelsize=14)
            plt.legend(title='Scaffold ID', fontsize=14)


    else:
        emsg = f"Unknown 'i_color_opt'={icolor_opt} in 'MultipleScf_Barplot' function. Admitted values: ['OneColor', 'Palette']"
        raise ValueError(emsg)
        
    # Set ticks on y axis
    new_y_ticks = range(1, max(iFreqDF["Frequency"])+1)  # Replace this list with your desired tick positions
    plt.yticks(new_y_ticks)

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

    plt.close()
# DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, SummaryPlotsDir, "OneColor", "red", "False", "--", 5, "ClusterSize_Distribution_by_Scaffolds.png", False)
# Change of visualization of zero frequency lines depending on the number of chromosomes (with visualization purpose only)
### if len(DF_Chr["ScfName"].unique()) < 3:
###    DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, SummaryPlotsDir, "OneColor", "blue", "False", "-", 5, "ClusterSize_Distribution_by_Scaffolds.svg", False)
### else:
###    DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, SummaryPlotsDir, "OneColor", "blue", "False", "--", 5, "ClusterSize_Distribution_by_Scaffolds.svg", False)


def ScfbyScf_Barplot(iDF_byChrFreqGenenum, iScaffoldIDcolname, iNumOfClusters, oDirname, icolor_opt, icolorOBJ, iVertical, iZeroLineStyle, iZeroLineSize, iFamName, iTwoFam_opt=False):
    DF_byChrFreqGenenum_gb = iDF_byChrFreqGenenum.groupby(iScaffoldIDcolname)

    for kScf, vDF in DF_byChrFreqGenenum_gb:
        # Out plot name and format
        if iFamName == "merged":
            Plotname = f"{iFamName}_" + kScf + ".svg"
        else:
            if iTwoFam_opt == False:
                Plotname = f"{iFamName}_fam_" + kScf + ".svg"
            else:
                Plotname = f"{iFamName}_2fam_" + kScf + ".svg"
        
        # Create a new dataframe, with 0 frequencies (with plotting purpose only)
        newDF = ModDf(vDF, iScaffoldIDcolname, iNumOfClusters)

        # Create a Barplot
        MultipleScf_Barplot(vDF, iScaffoldIDcolname, newDF, iNumOfClusters, oDirname, icolor_opt, icolorOBJ, iVertical, iZeroLineStyle, iZeroLineSize, iFamName, Plotname)
        
        
#### ScfbyScf_Barplot(DF_ChrFreqGenenum, ScaffoldIDcol, NumOfClusters, SinglePlotsDir, "OneColor", "pink", "AddVerticalLine", "-", 5)


# Save here the names of subdirs for each case
subdirs = ["merged"]

for iFamID, vDF in DF_1Fam.groupby("FamID"):
    subdirs.append(iFamID)
    
if len(DF_2Fam) != 0:
    temp = list(DF_2Fam["FamID"].unique())
    FamID = ".".join(temp)
    subdirs.append(FamID)

for idir in subdirs:
    if os.path.exists(f"{SummaryPlotsDir}/{idir}"):
        pass
    else:
        os.mkdir(f"{SummaryPlotsDir}/{idir}")

    if os.path.exists(f"{SinglePlotsDir}/{idir}"):
        pass
    else:
        os.mkdir(f"{SinglePlotsDir}/{idir}")



''' General plots considering the data from both families at the same time (as if they were one single family )'''

# Create a dataframe with the following information:
# Number of genes per cluster in each chr/scf
DF_Chr, ClusterSizes = GroupBy_ChrClusterID(DFcp, ScaffoldIDcol)

# Create a General Bar Plot (Frequency of cluster sizes)
# GeneralSummary_Plot(ClusterSizes, "lightgrey", "Global_ClusterSize_Distribution.png")
GeneralSummary_Plot(ClusterSizes, "lightgrey", "merged/Global_ClusterSize_Distribution_merged.svg", "merged")

# DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, SummaryPlotsDir, "OneColor", "red", "False", "--", 5, "ClusterSize_Distribution_by_Scaffolds.png", False)
# Change of visualization of zero frequency lines depending on the number of chromosomes (with visualization purpose only)
if len(DF_Chr["ScfName"].unique()) < 3:
    DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, f"{SummaryPlotsDir}/merged/", ColorOpt, ColorValueChoice, "False", "-", 5, "ClusterSize_Distribution_by_Scaffolds_merged.svg", "merged", False)
else:
    DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, f"{SummaryPlotsDir}/merged/", ColorOpt, ColorValueChoice, "False", "--", 5, "ClusterSize_Distribution_by_Scaffolds_merged.svg", "merged", False)

# Create Barplots
ScfbyScf_Barplot(DF_ChrFreqGenenum, ScaffoldIDcol, NumOfClusters, f"{SinglePlotsDir}/merged/", "OneColor", "lightgrey", "AddVerticalLine", "-", 5, "merged")



for iFamID, vDF in DF_1Fam.groupby("FamID"):
    
    # Create a dataframe with the following information:
    # Number of genes per cluster in each chr/scf
    DF_Chr, ClusterSizes = GroupBy_ChrClusterID(vDF, ScaffoldIDcol)

    # Create a General Bar Plot (Frequency of cluster sizes)
    GeneralSummary_Plot(ClusterSizes, "lightgrey", f"{iFamID}/{iFamID}_fam_Global_ClusterSize_Distribution.svg", iFamID)    
    
    # DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, SummaryPlotsDir, "OneColor", "red", "False", "--", 5, "ClusterSize_Distribution_by_Scaffolds.png", False)
    # Change of visualization of zero frequency lines depending on the number of chromosomes (with visualization purpose only)
    if len(DF_Chr["ScfName"].unique()) < 3:
        DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, f"{SummaryPlotsDir}/{iFamID}", ColorOpt, ColorValueChoice, "False", "-", 5, f"{iFamID}_fam_ClusterSize_Distribution_by_Scaffolds.svg", iFamID, False)
    
    else:
        DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, f"{SummaryPlotsDir}/{iFamID}", ColorOpt, ColorValueChoice, "False", "--", 5, f"{iFamID}_fam_ClusterSize_Distribution_by_Scaffolds.svg", iFamID, False)
    # Create Barplots
    ScfbyScf_Barplot(DF_ChrFreqGenenum, ScaffoldIDcol, NumOfClusters, f"{SinglePlotsDir}/{iFamID}", "OneColor", "lightgrey", "AddVerticalLine", "-", 5, iFamID)


if len(DF_2Fam) == 0:
    print("Warning! Two-family clusters dataframe has 0 length. No graphs will be printed")
else:
    temp = list(DF_2Fam["FamID"].unique())
    FamID = ".".join(temp)
    
    # Create a dataframe with the following information:
    # Number of genes per cluster in each chr/scf
    DF_Chr, ClusterSizes = GroupBy_ChrClusterID(DF_2Fam, ScaffoldIDcol)

    # Create a General Bar Plot (Frequency of cluster sizes)
    GeneralSummary_Plot(ClusterSizes, "lightgrey", f"{FamID}/{FamID}_2fam_Global_ClusterSize_Distribution.svg", FamID)    

    # DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, SummaryPlotsDir, "OneColor", "red", "False", "--", 5, "ClusterSize_Distribution_by_Scaffolds.png", False)
    # Change of visualization of zero frequency lines depending on the number of chromosomes (with visualization purpose only)
    if len(DF_Chr["ScfName"].unique()) < 3:
        DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, f"{SummaryPlotsDir}/{FamID}", ColorOpt, ColorValueChoice, "False", "-", 5, f"{FamID}_2fam_ClusterSize_Distribution_by_Scaffolds.svg", FamID, False)

    else:
        DF_ChrFreqGenenum, NumOfClusters = GeneralBar_Plot(DF_Chr, ScaffoldIDcol, f"{SummaryPlotsDir}/{FamID}", ColorOpt, ColorValueChoice, "False", "--", 5, f"{FamID}_2fam_ClusterSize_Distribution_by_Scaffolds.svg", FamID, False)
        
    # Create Barplots
    ScfbyScf_Barplot(DF_ChrFreqGenenum, ScaffoldIDcol, NumOfClusters, f"{SinglePlotsDir}/{FamID}", "OneColor", "lightgrey", "AddVerticalLine", "-", 5, FamID, True)
