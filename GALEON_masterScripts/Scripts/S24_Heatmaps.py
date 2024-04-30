import os, ast, sys, shutil
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
import itertools
from itertools import combinations
from mpl_toolkits.axes_grid1 import make_axes_locatable
from string import ascii_lowercase

plt.rcParams.update({'figure.max_open_warning': 0})
np.seterr(divide = 'ignore') # Ignore this warning: RuntimeWarning: divide by zero encountered in log

evo_mode= sys.argv[1]

if evo_mode == "enabled":
    g_threshold = sys.argv[2] # Ex: 100
    # g_threshold_shortID = float(g_threshold) / 100
    g_threshold_shortID = float(g_threshold) # update (22 Abril 2024)
    genefamily = sys.argv[3] # "GR_fam.gff3_matrices" 

    wkd_temp = sys.argv[4] # "../PhysicalDist_Matrices/"
    MergedDist_outdir = f"{os.path.dirname(wkd_temp)}/MergedDistances_Dataframes/"
    
    EvoMatricex_dir = sys.argv[5] # "Proteins"
    i_EvoMatrixPosition = sys.argv[6] # "Lower" / "Upper"
    
    evo_matrix_sep = sys.argv[7] # "tab"
    legend_scaleunits = sys.argv[8] # "default" or "auto"
    legend_evodist_decimals = int(sys.argv[9]) # "int"
    legend_physdist_decimals = int(sys.argv[10]) # "int"
    
    ColorMap_option = sys.argv[11] # "one" # "two"
    square_color_1 = sys.argv[12] # "black" # "black" 
    Colormap_1 = sys.argv[13] # "jet_r"  # "Heatmap Color (if 1 scale is used): ": HCmap_uniq,
    Colormap_2 = sys.argv[14] # "viridis_r"

    FamNum_opt = "WithinFamilies"


    scale_value = 1_000_000 # units: "1 Mb"

elif evo_mode == "disabled":
    g_threshold = sys.argv[2] # Ex: 100
    # g_threshold_shortID = float(g_threshold) / 100
    g_threshold_shortID = float(g_threshold) # update (22 Abril 2024)
    genefamily = sys.argv[3] # "GR_fam.gff3_matrices" 
    
    wkd_temp = sys.argv[4] # "../PhysicalDist_Matrices/"
    
    FamNum_opt = sys.argv[5]
    
    ColorMap_option = sys.argv[6] # "one" # "two"
    square_color_1 = sys.argv[7] # "black" # "black" # FAM1 cluster color
    square_color_2 = sys.argv[8] # "red" # "black"  # FAM2 cluster color
    square_color_3 = sys.argv[9] # "red" # "black" # FAM1+FAM2 cluster color 
    
    Colormap_1 = sys.argv[10] # "jet_r"  # "Heatmap Color (if 1 scale is used): ": HCmap_uniq,
    Colormap_2 = sys.argv[11] # "viridis_r" 

    scale_value = 1_000_000 # units: "1 Mb"

    legend_scaleunits = sys.argv[12] # "default" or "auto"
    legend_physdist_decimals = int(sys.argv[13]) # "int"

if ColorMap_option == "one":
    i_CMAP_list = [Colormap_1, Colormap_1]
elif ColorMap_option == "two":
    i_CMAP_list = [Colormap_1, Colormap_2]
else:
    msg = f"! Error: unknown 'ColorMap_option' {ColorMap_option}. Allowed values: 'one' or 'two'"
    raise ValueError(msg)
    
########## Input data ###############

def find_inputfiles(i_wkd, i_genefamily):
    ''' Function to detect dictionary files to be processed'''
    
    filelist = os.listdir(i_wkd)
    dictfile = ""
    # print(filelist, i_genefamily)

    for i in filelist:

        if i.endswith(".Dict.cluster.dict"):
            dictfile = i

            print(f"# ~~ Input dict file: {dictfile}")
            dictfile = f"{i_wkd}/{dictfile}"
        
        else:
            pass
            
    return dictfile


def discard_samefamily_clusters(i_dictfile, iD_clusters):
    odict_name = i_dictfile.replace("cluster.dict", "bw.cluster.dict")
    ofamlist = set()

    D_retained_info = {}

    for Scf, v1 in iD_clusters.items():
        for ilist in v1:
            genelist = ilist[-1]
            famlist = [ii.split("_")[0] for ii in genelist]
            famlist_set = set(famlist)
            ofamlist.update(famlist_set)

            if len(famlist_set) != 1:
                outinfo = ilist

                if Scf not in D_retained_info.keys():
                    D_retained_info[Scf] = []
                    D_retained_info[Scf].append(outinfo)
                else:
                    D_retained_info[Scf].append(outinfo)
                    
    # with open(odict_name, "w") as o1: # Sept 1 2023
        # o1.write(str(odict_name)) # Sept 1 2023

    return D_retained_info, list(ofamlist)

def import_dict(i_dictfile, i_filter=False):
    ''' Function for importing a dictionary from a text file'''
    
    with open(i_dictfile) as ftemp:
        
        xtemp = ftemp.read()
        odict = ast.literal_eval(xtemp)

        if i_filter == True:
            odict_filtered, ofamIDs = discard_samefamily_clusters(i_dictfile, odict)
            return odict, odict_filtered, ofamIDs

        elif i_filter == False:
            return odict

        else:
            raise ValueError()
    
    
def convert_sep_text(i_sep):
    if i_sep == "tab":
        i_sep = "\t"
        return i_sep
    
    elif i_sep == "comma":
        i_sep = ","
        return i_sep
    
    elif i_sep == "space":
        i_sep = " "
        return i_sep
    
    else:
        return i_sep
    
    
def cr_dir(i_dir, replace_opt=True): 
    if replace_opt == True:
        # Fx to create/overwrite a directory
        try:
            os.mkdir(i_dir)
                
        except FileExistsError:
            shutil.rmtree(i_dir)
            os.mkdir(i_dir)
    
    elif replace_opt == False:
        if os.path.exists(i_dir) == True:
            pass
        else:
            os.mkdir(i_dir)
        
class CustomTicker_default():
  
    def __call__(self, x, pos = None):
        
        if x in [np.nan, np.inf, -np.inf]:
            return 0
        
        else:
            orig_value = np.exp(x)
            # print("original value: ",orig_value)

            if orig_value < 1:
                output = round(orig_value, 2)
                return output
            else:
                temp = round(orig_value, 0)
                output = int(temp)
                return output

class CustomTicker_auto():
    def __init__(self, decimals=4):
        self.decimals = decimals
  
    def __call__(self, x, pos=None):
        if x in [np.nan, np.inf, -np.inf]:
            return 0
        else:
            orig_value = np.exp(x)
            output = round(orig_value, self.decimals)
            return output

class CustomTicker_special():
  
    def __call__(self, x, pos = None):
        
        if x in [np.nan, np.inf, -np.inf]:
            return 0
        
        else:
            orig_value = np.exp(x)
            # print("original value: ",orig_value)

            if orig_value < 1:
                output = round(orig_value, 4)
                return output
            else:
                temp = round(orig_value, 0)
                output = int(temp)
                return output

########## Evolutionary and Physical Heatmap Functions ###############
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def check_and_load_evo_input(i_famname, i_obj, i_EvoFieldSeparator):
    
    # Create a dict that associates each family to its Evomatrix
    i_famname = get_fam_name(i_famname)
    out_edist_dict = {}

    if os.path.isdir(i_obj) == True:

        tempdir = f"{os.getcwd()}/{i_obj}"
        temp_filelist = os.listdir(tempdir)
        evodist_tsvfile = []

        for i in temp_filelist:

            if i.startswith(i_famname) and i.endswith(".tsv"):
                evodist_tsvfile.append(i)
            else:
                pass

        if len(evodist_tsvfile) == 1:
            
            i_file = f"{tempdir}/{evodist_tsvfile[0]}"
            # out_edist_dict = import_evo_matrices(i_famname, i_file, i_EvoFieldSeparator)

            out_edist_dict[i_famname] = pd.read_csv(i_file, sep=i_EvoFieldSeparator, header=0, index_col=0)
            return out_edist_dict

        else:
            if len(evodist_tsvfile) == 0:
                msg = f"No evolutionary distance *tsv file has been found in this directory:\n- {tempdir}"
                raise ValueError(msg)

            elif len(evodist_tsvfile) > 1:
                msg = f"More than one *tsv file for the same family ('{i_famname}') has been found in this directory:\n- {tempdir}."
                raise ValueError(msg)
            else:
                raise ValueError("Unknown error!")
    
    elif os.path.isfile(i_obj) == True:
        i_file = i_obj
        
        # out_edist_dict = import_evo_matrices(i_famname, i_file, i_EvoFieldSeparator)
        out_edist_dict[i_famname] = pd.read_csv(i_file, sep=i_EvoFieldSeparator, header=0, index_col=0)
        return out_edist_dict
    
    else:
        msg = f"Input object: '{i_obj}' doesn't exist!"
        raise ValueError(msg)

def get_fam_name(i_text):
    # Get family name from phys matrix name
    
    fam_name = i_text.split(".")[0]
    return fam_name


def Merge_distances(i_Phys_Df, i_Evo_Df, i_scaleval, i_substitution, o_file):
    ''' Merge two distance matrices'''
    
    p_df_copy = i_Phys_Df.copy(deep=True)
    p_df_copy = p_df_copy.div(i_scaleval) # Scale the phys. dist. values

    Pcolnames = list(p_df_copy.columns)

    # Gene pairs combo
    Pgene_comb = list(combinations(Pcolnames, 2))
    
    # Out tabular format file: gene1, gene2, phys_dist, evo_dist (same data as in output dataframe)

    o_dictfile = o_file + ".distances.tsv"

    with open(o_dictfile, "w") as D_outfile:
        print("Gene1", "Gene2", "PhysicalDistance", "EvolutionaryDistance", sep="\t", end="\n", file=D_outfile)

        for i_gpair in Pgene_comb:
            try:
                Gene1, Gene2 = i_gpair[0], i_gpair[1]
                
                # Get values: Phys dist and Evo dist
                G1G2_physdist = i_Phys_Df.loc[Gene2,Gene1]
                G1G2_evodist = i_Evo_Df.loc[Gene1,Gene2]

                # Replace phys dist value by evo value
                p_df_copy.at[Gene2,Gene1] = G1G2_evodist
                
                print(Gene1, Gene2, G1G2_physdist, G1G2_evodist, sep="\t", end="\n", file=D_outfile)
                # print(Gene1, Gene2, G1G2_physdist, G1G2_evodist)
            except KeyError as e:
                sys.exit(f"! Key Error: {e}. Consider to check if your protein names match the gene names")

    if i_substitution == "Lower":

        o_file = o_file + ".EvoDist-Low_PhysDist-Up.matrix"

        # write to file
        p_df_copy.to_csv(o_file, sep="\t")

        return p_df_copy, o_file

    elif i_substitution == "Upper":

        o_file = o_file + ".EvoDist-Up_PhysDist-Low.matrix"

        # write to file
        p_df_copy.T.to_csv(o_file, sep="\t")

        return p_df_copy.T, o_file

    else:
        raise ValueError("# Unknown error in 'Merge_distance' function")

# Other functions
def create_alphanames():
    for size in itertools.count(1):
        for s in itertools.product(ascii_lowercase, repeat=size):
            yield "".join(s)

def get_alphanames(i_list):
    temp = []
    D_gene_alphabet_names = {}
    i_len = len(i_list)
    
    # get alphabet list
    for s in create_alphanames():
        temp.append(s)
        if len(temp) == i_len:
            break
    
    # join genename with alphabet
    for a,b in zip(i_list, temp):
        aname = f"{a} (# {b})"
        # D_gene_alphabet_names.append(aname)
        D_gene_alphabet_names[a] = aname
    
    temp = [f"# {ii}" for ii in temp]
    return D_gene_alphabet_names, temp




def ParseInputMatrix(i_matrixfile, i_matrixlocation):
    PlotOptionsHint = "Normal"

    Merged_df = pd.read_csv(f"{i_matrixlocation}/{i_matrixfile}", sep="\t", index_col=0) # Load df
        
    P_colnames = list(Merged_df.columns) # Get gene names (= column names)
    totalgenes = len(P_colnames)

    #::::# Graph parameteres
    # Out figure parameteres
    dpi_param = 200
    squaresize = 100 # pixels

    # Total number of genes
    n = totalgenes
    m = totalgenes

    if totalgenes < 15:
        n += 15
        m += 15
        PlotOptionsHint = "Special"
    else:
        pass

    heigh_width_ratio = n/m

    # set width of empty column used to stretch layout
    emptycol = 0.02
    # r = np.array([heigh_width_ratio, emptycol])

    # Figure height, width
    figwidth = m*squaresize/float(dpi_param)
    figheight = n*squaresize/float(dpi_param)

    # Borders
    left = 0.4; right = 0.6
    bottom = 0.4; top = 0.6

    # wspace (average relative space between subplots)
    wspace = 0.1

    #::::# Define Grid spaces
    Grid_plot = gridspec.GridSpec(1, 1)

    #::::# Plot adjustment
    plt.figure(figsize=(figwidth,figheight), dpi=dpi_param)

    #::::# Create a Plot
    EmptyCanvas = plt.subplot(Grid_plot[0, 0]) # Heatmap
    
    return Merged_df, P_colnames, EmptyCanvas, PlotOptionsHint


def Plot_and_Save_PhysDistOnly(i_matrixname, i_canvas, i_df, i_clusterdict, i_phys_threshold, i_autoUnitsScale, i_PhysScale_decimals,
                                o_directory, i_cmap, i_dfmaxvalue, i_pad=0.5, i_square_color_1="black", i_square_color_2="blue", i_square_color_3="red",
                                i_FamNum_opt="WithinFamilies", i_FilteredClusterDict=None, i_df_annot=False, dpi_param=200):
    

    # Create upper/lower semi-matrices
    mask1 = np.tril(np.ones_like(i_df, dtype=bool))
    mask2 = np.triu(np.ones_like(i_df, dtype=bool))
    
    #::::# Create additional axes for colorbars
    divider = make_axes_locatable(i_canvas)
    cAx2 = divider.append_axes("right", size = "2%", pad=i_pad)

    cAx2.set(xlabel=None)

    if i_autoUnitsScale == "default":
        if i_dfmaxvalue > 0.1:
            #::::# Create heatmaps
            smap44 = sns.heatmap(i_df, \
                    mask=mask1, \
                    annot = i_df_annot, \
                    square = True, \
                    robust = True, \
                    ax=i_canvas, cbar_ax=cAx2, \
                    cmap=i_cmap, \
                    cbar_kws={'ticks' : i_phys_threshold, 'format': CustomTicker_default()})

            smap45 = sns.heatmap(i_df, \
                    mask=mask2, \
                    annot = i_df_annot, \
                    robust = True, \
                    # square = True, \
                    ax=smap44, cbar=False, \
                    cmap=i_cmap)
        else:
            #::::# Create heatmaps
            smap44 = sns.heatmap(i_df, \
                    mask=mask1, \
                    annot = i_df_annot, \
                    square = True, \
                    robust = True, \
                    ax=i_canvas, cbar_ax=cAx2, \
                    cmap=i_cmap, \
                    cbar_kws={'format': CustomTicker_special()})

            smap45 = sns.heatmap(i_df, \
                    mask=mask2, \
                    robust = True, \
                    annot = i_df_annot, \
                    # square = True, \
                    ax=smap44, cbar=False, \
                    cmap=i_cmap)

    elif i_autoUnitsScale == "auto":
        #::::# Create heatmaps
        smap44 = sns.heatmap(i_df, \
                mask=mask1, \
                annot = i_df_annot, \
                square = True, \
                robust = True, \
                ax=i_canvas, cbar_ax=cAx2, \
                cmap=i_cmap, \
                cbar_kws={'format': CustomTicker_auto(decimals=i_PhysScale_decimals)})

        smap45 = sns.heatmap(i_df, \
                mask=mask2, \
                robust = True, \
                annot = i_df_annot, \
                ax=smap44, cbar=False, \
                cmap=i_cmap)
    else:
        raise ValueError("Unknown 'i_autoUnitsScale' value. Choose from ['default', 'auto']")
    
    # Add rectangles to the heatmap
    # millorar això
    i_matrixname = i_matrixname.replace(".merged.matrix", ".matrix")
    i_clusterinfo = i_clusterdict[i_matrixname]
    
    Pcolnames = list(i_df.columns)

    if len(i_clusterinfo) == 0:
        pass
    else:
        if i_FamNum_opt == "WithinFamilies":
            for i_entry in i_clusterinfo:
                # print(i_entry)
                c_shape = i_entry[2]
                c_keys = i_entry[3]
                c_keys_index = [Pcolnames.index(i) for i in c_keys]
                i_canvas.add_patch(Rectangle((c_keys_index[0],c_keys_index[0]), c_shape, c_shape, fill=False, linewidth=4, edgecolor=i_square_color_1))
        elif i_FamNum_opt == "BetweenFamilies":

            for i_entry in i_clusterinfo:
                c_shape = i_entry[2]
                c_keys = i_entry[3]
                c_keys_index = [Pcolnames.index(i) for i in c_keys]

                if i_matrixname in i_FilteredClusterDict.keys():
                    if i_entry in i_FilteredClusterDict[i_matrixname]:
                        i_canvas.add_patch(Rectangle((c_keys_index[0],c_keys_index[0]), c_shape, c_shape, fill=False, linewidth=4, edgecolor=i_square_color_3))
                    else:
                        famIDhint = GetNum_of_Fams(i_entry, "getfamname")
                        i_canvas.add_patch(Rectangle((c_keys_index[0],c_keys_index[0]), c_shape, c_shape, fill=False, linewidth=4, edgecolor=DictColor_bw2fam[famIDhint]))

                else:
                    famIDhint = GetNum_of_Fams(i_entry, "getfamname")
                    i_canvas.add_patch(Rectangle((c_keys_index[0],c_keys_index[0]), c_shape, c_shape, fill=False, linewidth=4, edgecolor=DictColor_bw2fam[famIDhint]))

        else:
            emsg = f"Unknown option '{i_FamNum_opt}'"
            raise ValueError(emsg)

    # Heatmap title
    i_canvas.set_title(f"{i_matrixname}")
    # Define labels for distance axes
    D_Pcolnames_alpha, alpha = get_alphanames(Pcolnames)
    Pcolnames_alpha = list(D_Pcolnames_alpha.values())

    if len(Pcolnames) <= 15:
        i_canvas.set_yticks(np.arange(len(Pcolnames)) + 0.5)
        i_canvas.set_yticklabels(Pcolnames_alpha, rotation=0, fontsize=10)

        i_canvas.set_xticks(np.arange(len(Pcolnames)) + 0.5)
        i_canvas.set_xticklabels(alpha, rotation = 90, fontsize=10)

        cAx2.tick_params(labelsize=10)

    elif len(Pcolnames) <= 30:
        i_canvas.set_yticks(np.arange(len(Pcolnames)) + 0.5)
        i_canvas.set_yticklabels(Pcolnames_alpha, rotation=0, fontsize=16)

        i_canvas.set_xticks(np.arange(len(Pcolnames)) + 0.5)
        i_canvas.set_xticklabels(alpha, rotation = 90, fontsize=16)

        cAx2.tick_params(labelsize=14)

    else:
        i_canvas.set_yticks(np.arange(len(Pcolnames)) + 0.5)
        i_canvas.set_yticklabels(Pcolnames_alpha, rotation=0, fontsize=20)

        i_canvas.set_xticks(np.arange(len(Pcolnames)) + 0.5)
        i_canvas.set_xticklabels(alpha, rotation = 90, fontsize=20)

        cAx2.tick_params(labelsize=20)

    # Export the figures
    Hfig33 = smap45.get_figure()
    Hfig33.figure.tight_layout(pad=2)
    Hfig33.savefig(f"{o_directory}/{i_matrixname}.PhysicalDist_Heatmap_{g_threshold_shortID}g.svg", dpi=dpi_param)
    Hfig33.savefig(f"{o_directory}/{i_matrixname}.PhysicalDist_Heatmap_{g_threshold_shortID}g.pdf", dpi=dpi_param)

    Hfig33.clf()


def Plot_and_Save_PhysEvoDist(i_matrixname, i_canvas, i_df1, i_df2, iPcolnames, iD_Pcolnames_alpha, i_clusterdict, 
                  i_phys_threshold, i_evo_threshold, i_autoUnitsScale, i_Scale1_decimals, i_Scale2_decimals,
                  i_cmap1, i_cmap2, o_directory, i_square_color_1="black",
                  i_df1_annot=False, i_df2_annot=False, dpi_param=300):
    
    # Create upper/lower semi-matrices
    mask1 = np.tril(np.ones_like(i_df1, dtype=bool))
    mask2 = np.triu(np.ones_like(i_df2, dtype=bool))
    
    #::::# Create additional axes for colorbars
    divider = make_axes_locatable(i_canvas)
    cAx2 = divider.append_axes("right", size = "2%", pad=0.5)
    cAx3 = divider.append_axes("bottom", size = "2%", pad=0.5)

    cAx2.set(xlabel=None)
    cAx3.set(xlabel=None)


    #::::# Create heatmaps
    if i_autoUnitsScale == "default":
        smap44 = sns.heatmap(i_df1, \
                mask=mask1, \
                annot = i_df1_annot, \
                robust = True, \
                square = True, \
                ax=i_canvas, cbar_ax=cAx2, \
                cmap=i_cmap1, \
                cbar_kws={'ticks' : i_phys_threshold, 'format': CustomTicker_default()})

        smap45 = sns.heatmap(i_df2, \
                mask=mask2, \
                annot = i_df2_annot, \
                robust = True, \
                # square = True, \
                ax=smap44, cbar_ax=cAx3, \
                cmap=i_cmap2, \
                cbar_kws={'ticks' : i_evo_threshold, 'format': CustomTicker_default(), \
                            'orientation' : 'horizontal'})

    elif i_autoUnitsScale == "auto": # update 30 abril
        smap44 = sns.heatmap(i_df1, \
                mask=mask1, \
                annot = i_df1_annot, \
                square = True, \
                robust = True, \
                ax=i_canvas, cbar_ax=cAx2, \
                cmap=i_cmap1, \
                cbar_kws={'format': CustomTicker_auto(decimals=i_Scale1_decimals)})

        smap45 = sns.heatmap(i_df2, \
                mask=mask2, \
                robust = True,
                annot = i_df2_annot, \
                ax=smap44, cbar_ax=cAx3, \
                cmap=i_cmap2, \
                cbar_kws={'format': CustomTicker_auto(decimals=i_Scale2_decimals), \
                            'orientation' : 'horizontal'})
    else:
        raise ValueError("Unknown 'i_autoUnitsScale' value. Choose from ['default', 'auto']")

    # Invert Evo dist legend axis
    cAx3.invert_xaxis() # 


    # Remove overlapping genes from alphanames dictionaries
    discarted_c_keys = [ki for ki in iD_Pcolnames_alpha.keys() if ki.startswith("GG#")]
    if len(discarted_c_keys) != 0:
        for k in discarted_c_keys:
            del iD_Pcolnames_alpha[k]

    # Add rectangles to the heatmap
    # millorar això
    i_matrixname = i_matrixname.replace(".merged.EvoDist-Low_PhysDist-Up.matrix", "")
    i_matrixname = i_matrixname.replace(".merged.matrix", ".matrix")
    i_clusterinfo = i_clusterdict[i_matrixname]
    if len(i_clusterinfo) == 0:
        pass
    else:
        for i_entry in i_clusterinfo:
            # Square shape
            c_shape = i_entry[2]

            # Overlapping genes will not be displayed in Phys. dist vs Evo. dist plot
            c_keys = [ci for ci in i_entry[3] if "GG#" not in ci]
            rm_c_keys = [ci for ci in i_entry[3] if "GG#" in ci] # remove keys, specifically from this cluster
            
            # Correct the shape size
            c_shape = i_entry[2] - len(rm_c_keys)
            

            if len(c_keys) > 1: # Add a shape only if the number of clustered genes is > 2 (After removing overlapping genes)
                c_keys_index = [iPcolnames.index(i) for i in c_keys]
                i_canvas.add_patch(Rectangle((c_keys_index[0],c_keys_index[0]), c_shape, c_shape, fill=False, linewidth=4, edgecolor=i_square_color_1))

    # Heatmap title
    i_canvas.set_title(f"{i_matrixname}")
    
    # Define labels for distance axes
    Pcolnames_alpha = list(iD_Pcolnames_alpha.values())
    alpha = ["(#" + a.split("(#")[1] for a in Pcolnames_alpha]

    # 
    i_canvas.set_yticks(np.arange(len(iPcolnames)) + 0.5)
    i_canvas.set_yticklabels(Pcolnames_alpha, rotation=0)

    i_canvas.set_xticks(np.arange(len(iPcolnames)) + 0.5)
    i_canvas.set_xticklabels(alpha, rotation = 90)

    # Export the figures
    Hfig33 = smap45.get_figure()
    Hfig33.figure.tight_layout(pad=2)

    Hfig33.savefig(f"{o_directory}/{i_matrixname}.PhysicalEvolutionaryDist_Heatmap_{g_threshold_shortID}g.svg", dpi=dpi_param)
    Hfig33.savefig(f"{o_directory}/{i_matrixname}.PhysicalEvolutionaryDist_Heatmap_{g_threshold_shortID}g.pdf", dpi=dpi_param)

    Hfig33.clf()


def GetNum_of_Fams(i_elem, i_opt, i_color1="green", i_color2="blue"):
    if i_opt == "assign":
        # Out dict
        D_out = {}

        # Get genelist
        genelist = []
        for k, vinfo in i_elem.items():
            for ii in vinfo:
                genelist += ii[-1]

        # Infer the number of gene families
        genefam_ids = [i.split("fam_")[0] for i in genelist]
        genefam_ids = sorted(list(set(genefam_ids)))

        if len(genefam_ids) == 2:
            # Out dict
            D_out[genefam_ids[0]+"fam_"] = i_color1
            D_out[genefam_ids[1]+"fam_"] = i_color2

            return D_out
        else:
            print(i_elem, genefam_ids)
            raise ValueError("Something is wrong with in 'BetweenFamilies' analysis")
    
    elif i_opt == "getfamname":
        clustered_geneslist = i_elem[-1]
        temp = list(set([i.split("fam_")[0] for i in clustered_geneslist]))
        famID = temp[0] + "fam_"

        if len(temp) == 1:
            return famID
        else:
            raise ValueError(temp)

###################################
''' Processing '''

### Define paths
# Path to physical matrices directory
PhysMatrices_dir = f"{wkd_temp}/{genefamily}"

# Path to cluster search results
ClustInfo_dir = f"{wkd_temp}/{genefamily}_{g_threshold_shortID}g"

# Import Cluster search results
dfile = find_inputfiles(ClustInfo_dir, genefamily)

# Import physical distance dictionary
if FamNum_opt == "BetweenFamilies":
    clusterdict, clusterdict_filtered, famIDs = import_dict(dfile, True)
    if len(clusterdict_filtered) != 0:
        DictColor_bw2fam = GetNum_of_Fams(clusterdict_filtered, "assign", square_color_1, square_color_2)
    else:
        DictColor_bw2fam = {}
        if len(famIDs) != 2:
            raise ValueError("Unknown error")
        else:
            for iFAM, ycolor in zip(famIDs, [square_color_1, square_color_2]):
                iFAM = iFAM + "_"
                DictColor_bw2fam[iFAM] = ycolor
        
elif FamNum_opt == "WithinFamilies":
    clusterdict = import_dict(dfile, False)

else:
    emsg = f"Error. in 'FamNum_opt' value: {FamNum_opt}. Admitted values: ['WithinFamilies', 'BetweenFamilies']"
    raise ValueError(emsg)

# Import all the Evo. Distance Matrices as pandas dataframes and save them in a dictionary
if evo_mode == "enabled":
    evo_matrix_sep = convert_sep_text(evo_matrix_sep) # convert "tab" input to "\t" separator. An
    print(genefamily, EvoMatricex_dir)
    D_Evodist_matrices = check_and_load_evo_input(genefamily, EvoMatricex_dir, evo_matrix_sep)

    cr_dir(MergedDist_outdir, False)

    tempname = genefamily.split("_")[0]
    MergedDist_out_mxdir = f"{MergedDist_outdir}/{tempname}_fam.merged.matrices"
    cr_dir(MergedDist_out_mxdir, False)
    
    MergedDist_out_plotdir = f"{MergedDist_outdir}/{tempname}_fam.plots_{g_threshold_shortID}g"
    cr_dir(MergedDist_out_plotdir, True)

    IntermediateFiles_out_mxdir = f"{MergedDist_outdir}/{tempname}_fam.IntermediateFiles"
    cr_dir(IntermediateFiles_out_mxdir, True)
else:
    pass


######### EXECUTION ############

#::::# Define threshold for both, physical and evolutinary distances. This is the default preset of the distances

# Define cutoffs for physical distances
a = list(np.arange(0,1,0.2))
b = list(range(1,12,2))
x = list(range(20,120,20))
y = list(range(0,1200,200))[1:]

z_phys = a + b + x + y

# Define cutoffs for evo distances
a = list(np.arange(0,0.2,0.04))
b = list(np.arange(0.2,1,0.2))
c = list(range(1,10,1))

# a = list(np.arange(0,1,0.2))
# b = list(range(1,10,2))

z_evo = a + b + c

# Convert phys/evo distances to 'log e' values
z_phys_ln = [np.log(i) for i in z_phys]
z_evo_ln = [np.log(i) for i in z_evo]


def get_max_dist(DF):
    DF = DF.set_index("Unnamed: 0")
    return DF.max().max()


# Create heatmaps based on physical distances
for scffile, clusterinfo in clusterdict.items():
    print(scffile)
    print("# ~~ Input matrix: ", scffile)

    # Merge Data

    if evo_mode == "enabled":
        Pdist_df, Pcolnames, A1, plothint = ParseInputMatrix(scffile, PhysMatrices_dir)
            
        totalgenes = len(Pcolnames)
        print(f"- Total genes: {totalgenes} in scf")
            
        ### Convert "merged matrix" units to log e base
            
        # Convert units
        Pdist_df = Pdist_df.replace([np.inf, -np.inf, "?"], np.nan)
        Pdist_df = Pdist_df.astype(float)
        newdf = Pdist_df/scale_value
        dfmax_value = newdf.max().max()

        newdf = np.log(newdf) # Log e base            
            
        ### Plot the physical distances dataframe
        if plothint == "Normal":
            Plot_and_Save_PhysDistOnly(scffile, A1, newdf, clusterdict, z_phys_ln, legend_scaleunits, legend_physdist_decimals, ClustInfo_dir, i_CMAP_list[0], dfmax_value, 0.5, square_color_1)
        elif plothint == "Special":
            Plot_and_Save_PhysDistOnly(scffile, A1, newdf, clusterdict, z_phys_ln, legend_scaleunits, legend_physdist_decimals, ClustInfo_dir, i_CMAP_list[0], dfmax_value, 0.1, square_color_1)

        ### Plot the physical distances + evolutionary distances
            
        associated_evo_mx = get_fam_name(scffile)
        evo_matrix = D_Evodist_matrices[associated_evo_mx]
        print(scffile, evo_matrix.max().max())
        outname_df = f"{MergedDist_out_mxdir}/{scffile}"
            
        ## - Merge distances
        # Prepare a dictionary with alpha names for each gene
        D_Pcolnames_alpha, alpha = get_alphanames(Pcolnames)
        
        # Drop columns
        Col2Drop = [i for i in Pcolnames if "GG#" in i] # Skip overlapping genes

        if len(Col2Drop) != 0:
            Pdist_df = Pdist_df.drop(Col2Drop)
            Pdist_df = Pdist_df.drop(columns=Col2Drop)
            # p_df_copy = p_df_copy.copy(deep=True)
            Pcolnames = list(Pdist_df.columns)
        else:
            pass

        mDF, mDFname = Merge_distances(Pdist_df, evo_matrix, scale_value, i_EvoMatrixPosition, outname_df)
        # Pdist_df, Pcolnames, A1, plothint = ParseInputMatrix(scffile, PhysMatrices_dir)
        #::::# Create a Plot
        Grid_plot = gridspec.GridSpec(1, 1)
        A1 = plt.subplot(Grid_plot[0, 0]) # Heatmap
                
        # Convert units
        mDF = mDF.replace([np.inf, -np.inf, "?"], np.nan)
        mDF = mDF.astype(float)
        newdf = np.log(mDF) # Log e base
    
        ## --create a copy of 'mDF'
        mDFcp = newdf.copy(deep=True) # This Df has log e base
            
        ### Plot the dataframe
        
        if i_EvoMatrixPosition == "Lower":
            Plot_and_Save_PhysEvoDist(scffile, A1, newdf, mDFcp, Pcolnames, D_Pcolnames_alpha, clusterdict, 
                            z_phys_ln, z_evo_ln, legend_scaleunits, legend_physdist_decimals, legend_evodist_decimals,
                            i_CMAP_list[0], i_CMAP_list[1], MergedDist_out_plotdir, square_color_1)
        elif i_EvoMatrixPosition == "Upper":
            print(scffile, newdf.shape, mDFcp.shape)
            Plot_and_Save_PhysEvoDist(scffile, A1, newdf, mDFcp, Pcolnames, D_Pcolnames_alpha, clusterdict, 
                                z_evo_ln, z_phys_ln, legend_scaleunits, legend_evodist_decimals, legend_physdist_decimals, 
                                i_CMAP_list[0], i_CMAP_list[1], MergedDist_out_plotdir, square_color_1)

            
        
    elif evo_mode == "disabled":
        Pdist_df, Pcolnames, A1, plothint = ParseInputMatrix(scffile, PhysMatrices_dir)

        totalgenes = len(Pcolnames)
        print(f"- Total genes: {totalgenes} in scf")
            
        ### Convert "merged matrix" units to log e base
            
        # Convert units
        Pdist_df = Pdist_df.replace([np.inf, -np.inf, "?"], np.nan)
        Pdist_df = Pdist_df.astype(float)
        newdf = Pdist_df/scale_value
        dfmax_value = newdf.max().max()

        
        newdf = np.log(newdf) # Log e base            
            
        ### Plot the dataframe
        if FamNum_opt == "BetweenFamilies":
            if plothint == "Normal":
                Plot_and_Save_PhysDistOnly(scffile, A1, newdf, clusterdict, z_phys_ln, legend_scaleunits, legend_physdist_decimals, ClustInfo_dir, i_CMAP_list[0], dfmax_value, 0.5, square_color_1, square_color_2, square_color_3, FamNum_opt, clusterdict_filtered)
            elif plothint == "Special":
                Plot_and_Save_PhysDistOnly(scffile, A1, newdf, clusterdict, z_phys_ln, legend_scaleunits, legend_physdist_decimals, ClustInfo_dir, i_CMAP_list[0], dfmax_value, 0.1, square_color_1, square_color_2, square_color_3, FamNum_opt, clusterdict_filtered)
        
        elif FamNum_opt == "WithinFamilies":
            if plothint == "Normal":
                Plot_and_Save_PhysDistOnly(scffile, A1, newdf, clusterdict, z_phys_ln, legend_scaleunits, legend_physdist_decimals, ClustInfo_dir, i_CMAP_list[0], dfmax_value, 0.5, square_color_1, square_color_2, square_color_3, FamNum_opt, None)
            elif plothint == "Special":
                Plot_and_Save_PhysDistOnly(scffile, A1, newdf, clusterdict, z_phys_ln, legend_scaleunits, legend_physdist_decimals, ClustInfo_dir, i_CMAP_list[0], dfmax_value, 0.1, square_color_1, square_color_2, square_color_3, FamNum_opt, None)

    else:
        emsg = f"Unknown 'evo_mode' parameter: {evo_mode}"
        raise ValueError(emsg)
            
        

