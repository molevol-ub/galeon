import os, ast, gc, sys
import pandas as pd
from itertools import combinations
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
import matplotlib.ticker as mticker
from matplotlib.ticker import LogFormatterSciNotation
from matplotlib.ticker import MaxNLocator


''' Input Parameteres '''
# input g value
g_threshold = sys.argv[1] # 100 # arg 4
# g_threshold_shortID = float(g_threshold) / 100
g_threshold_shortID = float(g_threshold) # update (22 Abril 2024)


# Family physical matrices dirname
FamMatrix_dir = sys.argv[2] # "GR_fam.gff3_matrices" # arg 2
FamShortID = FamMatrix_dir.split("_")[0]

# Path to physical matrices
Pmatrices_dir = sys.argv[3] # "../PhysicalDist_Matrices" # arg 1

# Path to merged distances data
MergedDist_dir = sys.argv[4] # "../MergedDistances_Dataframes" # arg 3
MergedDist_dir = f"{os.path.dirname(Pmatrices_dir)}/{MergedDist_dir}" # "../MergedDistances_Dataframes" # arg 3
FamMergedDistData_dir = f"{FamShortID}_fam.merged.matrices" # Dir with files containing both distances

print(g_threshold, FamMatrix_dir, Pmatrices_dir, MergedDist_dir)

####
FamCase_dir_1 = FamMergedDistData_dir
FamCase_dir_2 = f"{FamMatrix_dir}_{g_threshold_shortID}g"

# Path to save the plots
FamPlot_dir = f"{FamShortID}_fam.plots_{g_threshold_shortID}g"
FamPlot_dir = f"{MergedDist_dir}/{FamPlot_dir}" # Scatterplots are going to be saved here

# Path to save intermediate files
IntermediateFiles_dir = f"{MergedDist_dir}/{FamShortID}_fam.IntermediateFiles" # Scatterplots are going to be saved here


# List all *distances.tsv files
DistDir = f"{MergedDist_dir}/{FamCase_dir_1}"
tsv_files = [i for i in os.listdir(DistDir) if i.endswith(".distances.tsv")]

# List all the *dict files
DictDir = f"{Pmatrices_dir}/{FamCase_dir_2}"
dict_files = [i for i in os.listdir(DictDir) if i.endswith(".dict")]
print(dict_files)

user_outformat="svg"



''' Load the dictionary containing cluster information '''

# Out figure name
dfile = f"{DictDir}/{dict_files[0]}"

# Import and save the dict with cluster information
D_Cluster_Info = {}

## - Load the dictionary
with open(dfile) as d:
    dfile_loaded = d.read()
    D_clust = ast.literal_eval(dfile_loaded)
    
# Create a new one with each gene pair set as a 'key' and each cluster ID set as 'value'
for matrixID, cinfo in D_clust.items():
    ScfID = matrixID.split("temp_matrices.")[-1].split(".matrix")[0]
    D_Cluster_Info[ScfID] = {}

    for ClusterCase in cinfo:
        ClusterID, ClusteredGenes = ClusterCase[0], set(ClusterCase[3])
        
        # Every possible combination of clustered genes
        Pgene_comb = list(combinations(ClusteredGenes, 2))
        
        
        for k in Pgene_comb:
            g1, g2 = k[0], k[1]
            hintID = f"{g1}-@-{g2}" # Label
            hintID_r = f"{g2}-@-{g1}" # Rv label
            
            # Add info to the dict
            D_Cluster_Info[ScfID][hintID] = ClusterID 
            D_Cluster_Info[ScfID][hintID_r] = ClusterID

# List to store all the dataframes
df_list = []

# TSV file
for f in tsv_files:
    # Path to input tsv file
    path_to_file = f"{DistDir}/{f}"
    
    temp = f.split("_matrices.")[1]
    ScfID = temp.split(".matrix.")[0]
    
    # Load file as dataframe
    out_headers = []
    out_list = []
    
    with open(path_to_file, "r") as f1:
        for line in f1:
            line = line.strip().split("\t")
            
            if "Dist" in line[2]:
                newheader = line
                newheader.append("ScfID")
                newheader.append("Hint")
                newheader.append("ClusterID")
                out_headers = newheader
            
            else:
                g1, g2 = line[0], line[1]
                hintID = '-@-'.join([g1, g2])

                newline = line
                newline.append(ScfID)
                newline.append(hintID)
                
                try:
                    clustID = D_Cluster_Info[ScfID][hintID]
                    newline.append(clustID)

                except KeyError:
                    clustID = None
                    newline.append("NC")
            
                out_list.append(newline)
                
                

    df = pd.DataFrame(out_list, index=None, columns=out_headers)
    df_list.append(df)
    
# Concatenate the loaded matrices
df_glob = pd.concat(df_list)
del df_list
gc.collect() # Free memory (useful when there is a huge amount of input data)


class LogFormatter_CustomTicker(LogFormatterSciNotation):
  
    def __call__(self, x, pos = None):
  
        if x not in [0.1, 1, 10]:
            return LogFormatterSciNotation.__call__(self, x, pos = None)
  
        else:
            # return "{x:g}".format(x = x)
            if x == 1:
                return 1
            elif x == 10:
                return 10
            else:
                return x

def Scf_Scatter_Plot(i_df, display_legend_opt):
    # Crete a plot
    # fig, A1 = plt.subplots()
    A1 = plt.gca()
    
    #print(type(A1))
    # Define x ticks
    # maximum allowed
    x = range(1,11)
    x = [10**i for i in x]
    
    # convert Distances columns  to numeric type
    i_df["PhysicalDistance"] = i_df["PhysicalDistance"].apply(pd.to_numeric)
    i_df.loc[i_df["EvolutionaryDistance"] == "?", "EvolutionaryDistance"] = None
    i_df["EvolutionaryDistance"] = i_df["EvolutionaryDistance"].apply(pd.to_numeric)
        
    # find the max physical distance from input df
    max_x = np.max(i_df["PhysicalDistance"])
    max_x_round = round(max_x, -3)
    
    if max_x_round == 0:
        max_x_round = max_x
        x_filtered = [xi for xi in x if xi <= max_x_round] # filter thresholds by discarding all those below a maximum
    else:
        x_filtered = [xi for xi in x if xi <= max_x_round] # filter thresholds by discarding all those below a maximum
        
    x_log = [np.log10(i) for i in x]
    
    # Plot
    Cids_unique = [i for i in list(i_df["ClusterID"].unique())]
    
    # ClusterID column mod
    i_df["Cluster_hint"] = i_df.loc[:, ["ClusterID"]]
    D_Cids = {i : i.replace("Cluster", "") for i in Cids_unique}
    D_Cids.update({"NC" : "-1"})
    # i_df["Cluster_hint"].replace(D_Cids, inplace=True)
    i_df["Cluster_hint"] = i_df["Cluster_hint"].replace(D_Cids)
    
    i_df["Cluster_hint"] = i_df["Cluster_hint"].apply(pd.to_numeric)
    i_df.sort_values(by=["Cluster_hint"],  inplace=True)
    
    mypalette = dict(zip(Cids_unique, sns.color_palette("gist_rainbow", n_colors=len(Cids_unique))))
    mypalette.update({'NC': (0.73, 0.73, 0.73)}) # Set NC color to grey

    sns_plot = sns.scatterplot(x="PhysicalDistance", y="EvolutionaryDistance", \
                               data=i_df, hue="ClusterID", ax=A1, \
                               alpha=0.9, s = 12, palette=mypalette, \
                               clip_on=True)

    # sns_plot = sns.scatterplot(x="phys", y="evo", data=i_df, hue="clusters", palette="rainbow", ax=A1, clip_on=True)
    A1.set_xlabel('Physical distances (bp)')
#     A1.set_ylabel('Evolutionary distances (aa. subst/position)')
    A1.set_ylabel('Evolutionary distances')
    
    # Display/Hide legend

    if display_legend_opt == True:

        leg = A1.legend()
        sns.move_legend(A1, "upper left", bbox_to_anchor=(1,1.02), fontsize = "medium")
        
    elif display_legend_opt == False:
        A1.get_legend().remove()
        
    # A1.set_xlim(1) # Definir on comença l'etiquetatge de l'eix X
    A1.set_ylim(-0.05) # Definir on comença l'etiquetatge de l'eix X
    A1.set_xscale("log", base=10) # Passar a escala log10
    # A1.set_xticks(x_filtered) # Definir els ticks
    
    #A1.xaxis.set_major_locator(MaxNLocator(nbins="auto"))
    A1.xaxis.set_major_formatter(LogFormatter_CustomTicker())
    
    print(x_filtered)
    return A1, i_df

# Plot each Scaffold separately
df_list_updated = []

for matrixID, cinfo in D_clust.items():
    
    temp = matrixID.split(".matrix")
    # scatter_plot_name = FamPlot_dir + "/" + ".".join(temp) + f".PhysicalEvolutionaryDist_Scatterplot_{g_threshold_shortID}g.svg"
    # scatter_plot_name2 = FamPlot_dir + "/" + ".".join(temp) + f".PhysicalEvolutionaryDist_Scatterplot_{g_threshold_shortID}g.pdf"

    ScfID = temp[0].split("temp_matrices.")[-1]

    scatter_plot_name = f'{FamPlot_dir}/{temp[0]}.matrix.PhysicalEvolutionaryDist_Scatterplot_{g_threshold_shortID}g.svg'
    scatter_plot_name2 = f'{FamPlot_dir}/{temp[0]}.matrix.PhysicalEvolutionaryDist_Scatterplot_{g_threshold_shortID}g.pdf'

    # Create a scatter plot for each scaffold
    t_df = df_glob.loc[df_glob["ScfID"] == ScfID].copy() # Select scaffold specific data

    print(t_df)
    if len(t_df) == 0:
        emsg = f"Empty dataframe for this scaffold: {ScfID}"
        sys.exit(emsg)

    t_plot, t_df_updated = Scf_Scatter_Plot(t_df, False) # Plot the graph
    df_list_updated.append(t_df_updated)
    
    # Export figure
    Sfig = t_plot.get_figure()
    Sfig.savefig(scatter_plot_name, dpi = 100, bbox_inches="tight", facecolor="white")
    Sfig.savefig(scatter_plot_name2, dpi = 100, bbox_inches="tight", facecolor="white")
    Sfig.clf()
    
''' Prepare the data to create a Global Plot with all the scaffolds '''

# Concatenate the loaded matrices
df_glob_updated = pd.concat(df_list_updated)
del df_list_updated
gc.collect()


# In[11]:


# Create a new column. It's a technical column with visualization purpose for NC (not clustered genes)
df_glob_updated["Scf_hint"] = df_glob_updated.loc[:, "ScfID"]
df_glob_updated.head()

# Update Scf_hint column based on Cluster_hint column
df_glob_updated.loc[df_glob_updated.Cluster_hint == -1, "Scf_hint"] = "NC"

# Save the dataframe 
OutFilename = f"{IntermediateFiles_dir}/{FamShortID}_Family_PhysDist_and_EvoDist_Raw_data4ScatterPlots.tsv"
print(OutFilename)
df_glob_updated.to_csv(OutFilename,sep="\t",index=None)


# In[12]:


def Glob_Scatter_Plot(i_df, display_legend_opt):
    # Crete a plot
    # fig, A1 = plt.subplots()
    A1 = plt.gca()
    
    # Define x ticks
    # maximum allowed
    x = range(1,11)
    x = [10**i for i in x]
    
    # convert Distances columns  to numeric type
    i_df["PhysicalDistance"] = i_df["PhysicalDistance"].apply(pd.to_numeric)
    i_df.loc[i_df["EvolutionaryDistance"] == "?", "EvolutionaryDistance"] = None
    i_df["EvolutionaryDistance"] = i_df["EvolutionaryDistance"].apply(pd.to_numeric)
        
    # find the max physical distance from input df
    max_x = np.max(i_df["PhysicalDistance"])
    max_x_round = round(max_x, -3)
    
    if max_x_round == 0:
        max_x_round = max_x
        x_filtered = [xi for xi in x if xi <= max_x_round] # filter thresholds by discarding all those below a maximum
    else:
        x_filtered = [xi for xi in x if xi <= max_x_round] # filter thresholds by discarding all those below a maximum
        
    x_log = [np.log10(i) for i in x]
    
    # Plot
    Cids_unique = [i for i in list(i_df["Scf_hint"].unique())]
    
    # Define palette of colors
    mypalette = dict(zip(Cids_unique, sns.color_palette("gist_rainbow", n_colors=len(Cids_unique))))
    mypalette.update({'NC': "#BABABA"}) # Set NC color to grey
    
    print(mypalette)
    sns_plot = sns.scatterplot(x="PhysicalDistance", y="EvolutionaryDistance", \
                               data=i_df, hue="Scf_hint", ax=A1, \
                               alpha=1, s = 25, palette=mypalette, \
                               clip_on=True)
    # sns_plot = sns.scatterplot(x="phys", y="evo", data=i_df, hue="clusters", palette="rainbow", ax=A1, clip_on=True)
    A1.set_xlabel('Physical distances (bp)')
    # A1.set_ylabel('Evolutionary distances (aa. subst/position)')
    A1.set_ylabel('Evolutionary distances')
    
    # Display/Hide legend

    if display_legend_opt == True:

        leg = A1.legend()
        sns.move_legend(A1, "upper left", bbox_to_anchor=(1,1.02), fontsize = "medium")
        
    elif display_legend_opt == False:
        A1.get_legend().remove()
        
    # A1.set_xlim(1) # Definir on comença l'etiquetatge de l'eix X
    A1.set_ylim(-0.05) # Definir on comença l'etiquetatge de l'eix X
    A1.set_xscale("log", base=10) # Passar a escala log10
    A1.set_ylim(0, 6)
    
    # A1.xaxis.set_major_locator(MaxNLocator(nbins="auto"))
    A1.xaxis.set_major_formatter(LogFormatter_CustomTicker())
    
    return A1

# Sort the dataframe by "Cluster_hint" column to set "NC" at the end of the legend
df_glob_updated.sort_values(by=["Cluster_hint"],  inplace=True)

glob_plot = Glob_Scatter_Plot(df_glob_updated, False)

# Export figure
dfile_basename_1 = os.path.basename(dfile).split(".")[0] + f".GlobScatterPlot_{g_threshold_shortID}g.svg"
dfile_basename_1 = f"{MergedDist_dir}/{dfile_basename_1}"
dfile_basename_2 = os.path.basename(dfile).split(".")[0] + f".GlobScatterPlot_{g_threshold_shortID}g.pdf"
dfile_basename_2 = f"{MergedDist_dir}/{dfile_basename_2}"

# Export dataframe
dfile_basename_3 = os.path.basename(dfile).split(".")[0] + f".GlobScatterPlot_{g_threshold_shortID}g_RawData.tsv"
dfile_basename_3 = f"{MergedDist_dir}/{dfile_basename_3}"
df_glob_updated.to_csv(dfile_basename_3, sep="\t", index=None)

Sfig = glob_plot.get_figure()
Sfig.savefig(dfile_basename_1, dpi = 100, bbox_inches="tight", facecolor="white")
Sfig.savefig(dfile_basename_2, dpi = 100, bbox_inches="tight", facecolor="white")
Sfig.clf()
