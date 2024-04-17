import ast, math, shutil, os, sys, argparse
from collections import Counter
import pandas as pd
from itertools import combinations
from scipy.stats import mannwhitneyu


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Final Report', formatter_class=argparse.RawTextHelpFormatter) #ArgumentDefaultsHelpFormatter)
    parser._optionals.title = "Input arguments"
    
    parser.add_argument("-clust", "--ClusterfinderResults",
        type = str,
        required = True,
        help = "(REQUIRED). GALEON clusterfinder results directory.")
    
    parser.add_argument("-coords", "--CoordsDir",
        type = str,
        required = True,
        help = "(REQUIRED) Input coordinates directory.")

    parser.add_argument("-prot", "--ProtDir",
        type = str,
        required = True,
        help = "Input protein directory")

    parser.add_argument("-outdir", "--OutDirectory",
        type = str,
        default = "MannWhitney_StatisticsResults",
        help = "Output directory.")
        
    # Parsing arguments
    args = parser.parse_args()


    # Get configvalue
    config = vars(args)
    print(config)


# In[ ]:


# Inputs
clusterfinder_dir = config["ClusterfinderResults"]
coords_dir = config["CoordsDir"]
protdir = config["ProtDir"]

# Outdir
MWoutdir = config["OutDirectory"]


# In[16]:


# # Inputs
# clusterfinder_dir = "clusterfinder_Results_Directory/"
# coords_dir = "GFFs/"
# protdir = "Proteins/"

# # Outdir
# MWoutdir = "MannWhitney_StatisticsResults/"


# In[17]:


# Input files
D_INPUTfiles = {}

for r, d, l in os.walk(f"{clusterfinder_dir}/PhysicalDist_Matrices/"):
    if len(l) != 0:
        dictlist = [i for i in l if i.endswith("cluster.dict")]

        if len(dictlist) == 0:
            pass
        
        elif len(dictlist) == 1:
            D_path = f"{r}/{dictlist[0]}"
            FAM_name = dictlist[0].split(".")[0]
            
            D_INPUTfiles[FAM_name] = {"BEDfile" : "", "EvoDistfile" : "", "ClusterDictfile" : []}
        else:
            raise ValueError("Something is wrong, multiple dictionaries were found")
            
for r, d, l in os.walk(f"{clusterfinder_dir}/PhysicalDist_Matrices/"):
    if len(l) != 0:
        dictlist = [i for i in l if i.endswith("cluster.dict")]

        if len(dictlist) == 0:
            pass
        
        elif len(dictlist) == 1:
            D_path = f"{r}/{dictlist[0]}"
            FAM_name = dictlist[0].split(".")[0]
            
            D_INPUTfiles[FAM_name]["ClusterDictfile"].append(D_path)
        else:
            raise ValueError("Something is wrong, multiple dictionaries were found")
            
# Search GFF(bed) files
I_GFF_list = [f"{coords_dir}/{i}" for i in os.listdir(f"{coords_dir}/") if i.endswith("collapsed.temp")]

# Search overlapping genes
I_Ovlap_list = [f"{coords_dir}/{i}" for i in os.listdir(f"{coords_dir}/") if i.endswith("OverlappingGenes.dict")]

# Search distances matrices
I_distmatrices_list = [f"{protdir}/{i}" for i in os.listdir(f'{protdir}/') if i.endswith('distancematrix.tsv')]
print(I_distmatrices_list)
# Assign input files
for kFAM, v in D_INPUTfiles.items():
    print(kFAM)
    for gfffile in I_GFF_list:
        if kFAM in gfffile:
            v["BEDfile"] = gfffile
            
    for evofile in I_distmatrices_list:
        if kFAM in evofile:
            v["EvoDistfile"] = evofile
            
    for ovlapfile in I_Ovlap_list:
        if kFAM in ovlapfile:
            v["OverlappingGenesFile"] = ovlapfile


# In[ ]:





# In[18]:


if os.path.exists(MWoutdir) == False:
    os.mkdir(MWoutdir)
else:
    shutil.rmtree(MWoutdir)
    os.mkdir(MWoutdir)


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


def Dck_fx(igene_num, iEvoDist_list):
    
    ifactor = 2/(igene_num*(igene_num - 1)) # 2/(n*(n-1))
    isummatory = sum(iEvoDist_list)

    Dck = ifactor*isummatory
    return Dck

def DT_fx(igene_num, iEvoDist_list):
    
    ifactor = 2/(igene_num*(igene_num - 1)) # 2/(n*(n-1))
    isummatory = sum(iEvoDist_list)

    DT = ifactor*isummatory
    return DT

def DT2_stat_glob(iDF_EvoDist):
    # Get combinations of all genes of a given family
    ALLgene_combo = list(combinations(iDF_EvoDist.index, 2))

    # Get evol distances for each gene pair
    ALLevodists = [iDF_EvoDist.loc[j[0],j[1]] for j in ALLgene_combo]
    
    # Get DT using over all the genome (considering all the scaffolds as one)
    DT2_glob = DT_fx(len(iDF_EvoDist.index), ALLevodists)
    
    return DT2_glob

def DC_fx(iclust_num, iDckValues_list):
    
    ifactor = 1/iclust_num
    isummatory = sum(iDckValues_list)
    
    DC = ifactor*isummatory
    
    return DC

def Get_Dckvalues(o_key, iClustInfo, iDF_EvoDist, o_dict):
    for iclust in iClustInfo:
        # Get 
        # - cluster id
        # - number of clustered genes
        # - list of clustered genes
        Cname, Cgenenum, Cgenelist = iclust[0], iclust[2], iclust[3]

        # Discard overlapping genes (starting with GG#...)
        discarted_genes = [i for i in Cgenelist if i.startswith("GG#") == True]

        if len(discarted_genes) != 0:
            Cgenelist = [i for i in Cgenelist if i not in discarted_genes]
            Cgenenum = Cgenenum - len(discarted_genes)

            for i in discarted_genes:
                print(f"Overlapping genes labeled as: '{i}', are not considered in evolutionary statistics computation")

        # Generate gene combinations
        Cgenecombo = list(combinations(Cgenelist, 2))

        # Get evol distances for each gene pair
        Cevodists = [iDF_EvoDist.loc[j[0],j[1]] for j in Cgenecombo]

        # Compute Dck for each cluster
        Dck_value = Dck_fx(Cgenenum, Cevodists)
        # print("Dck = ", Dck_value)
        
        # Export Dck_value
        o_dict[o_key]["Dck"][Cname] = Dck_value

def sigf_status(sgf_value):
    s_status = ""

    if sgf_value >= 0.05:
        s_status = "NS"
    elif sgf_value < 0.05 and sgf_value >= 0.01:
        s_status = "*"
    elif sgf_value < 0.01 and sgf_value >= 0.005:
        s_status = "**"
    elif sgf_value < 0.005:
        s_status = "***"

    return s_status

def other_statistics(n1, n2, U_statistic):

    # Compute: 'mw', 'variance' and 'std. deviation'
    mw_value =  (n1*n2)/2 # (n1*n2)/2
    variance_value = (n1*n2*(n1+n2+1))/12
    stdev_value = math.sqrt(variance_value)

    # Compute 'z score'

    z_score = (abs(U_statistic - mw_value)-0.5)/stdev_value
    
    return z_score, mw_value, variance_value, stdev_value

def get_gvalueID(i_dictname):
    outname = os.path.basename(i_dictname).split("matrices_")[-1].split(".Dict")[0]
    return outname

def Run_MannWhiteneyU(x,y):
    U1_Statistic, pvalue = mannwhitneyu(x, y, alternative="two-sided")
    U2_Statistic = len(x)*len(y) - U1_Statistic

    z_score, mw_value, variance_value, stdev_value = other_statistics(len(x), len(y), U1_Statistic)
    
    return [U1_Statistic, U2_Statistic, pvalue, z_score, mw_value, variance_value, stdev_value]


# In[84]:


def add_stats(i_dict, i_DF, i_statname):
    D_temp = {Scaffold.split("_matrices.")[1].split(".matrix")[0] : v[i_statname] for Scaffold, v in i_dict.items()}
    i_DF[i_statname] = i_DF["ScaffoldID"].map(D_temp)

def export_Dck(i_dict, o_file):
    outlist = []
    for Scaffold, vdict in i_dict.items():
        for k, v in vdict.items():
            if isinstance(v, dict):
                for kk, vv in v.items():
                    FamID= Scaffold.split("_fam.")[0] + "_fam"
                    ScaffoldID= Scaffold.split("_matrices.")[1].split(".matrix")[0]
                    outlist.append([FamID, ScaffoldID, k, kk, vv])

    DF = pd.DataFrame.from_records(outlist)
    # DF.columns = ["Family", "ScaffoldID", "Statistic", "ClusterID", "Statistic value"]
    DF.columns = ["FamID", "ScaffoldID", "Statistic", "ClusterID", "Statistic value"]
    
    DF.to_csv(o_file, sep="\t", index=None)

    
def GetStatistics_and_RunMWtest(i_FAMname, i_bedfile, i_evodistfile, i_clusterfile, i_gvhint, i_overlap_dictfile=None):

    ''' Load GFF collapsed bed file '''
    # The genes are already sorted by the start coordinate

    DF_FamBED = pd.read_csv(i_bedfile, sep="\t", header=None)
    # DF_FamBED.columns = ["ScaffoldID", "GeneStart", "GeneEnd", "Attribute", "bedclustID"]
    DF_FamBED.columns = ["ScaffoldID", "GeneStart", "GeneEnd", "Attribute", "bedclustID"]

    ''' Load Evolutionary Distance Matrix '''

    DF_EvoDist = pd.read_csv(i_evodistfile, sep="\t", index_col=0)

    ''' Load Cluster dict '''

    CLUSTERdictfile = i_clusterfile

    D_clusters = CheckClusterData(CLUSTERdictfile)
    
    if D_clusters == None:
        return "Empty dictionary"
    
    ''' Exclude overlapping entries (if any) '''
    # Load the information about overlapping genes
    if i_overlap_dictfile != None:
        D_ovlap = load_dict(i_overlap_dictfile)
        Genes2exclude = [] # Save the bed clustered IDs of overlapping genes
        BedClusterNames = [] # Save the IDs of overlapping genes
        
        for k,v in D_ovlap.items():
            BedClusterNames.append(k)
            for gene_entry in v:
                geneID = gene_entry[-2]
                Genes2exclude.append(geneID)

        print(f"- Warning: Overlapping genes won't be considered: {[Genes2exclude]}")

        # Exclude overlapping genes from the evo matrix
        DF_EvoDist.drop(columns=Genes2exclude, index=Genes2exclude, inplace=True)

        # Exclude all the entries with BedClusteredNames
        temp = DF_FamBED[DF_FamBED["Attribute"].isin(BedClusterNames)==False].copy(deep=True)
        DF_FamBED = temp

    # Empty dict to fill with the genelist IDs of each scaffold
    D_ScfFamGenes = {}

    # Group by scaffold and get the gene names

    DF_FamBED_gb = DF_FamBED.groupby("ScaffoldID")
    for Scf, Genes in DF_FamBED_gb:
        D_ScfFamGenes[Scf] = []

        for a in list(Genes["Attribute"]):
            geneID = a.split(";")[0].replace("ID=","")
            D_ScfFamGenes[Scf].append(geneID)

    ''' For each scaffold, compute Dck values for each input cluster list '''

    # Empty dict to store Dck values
    D_clust_Dck = {}

    for k,v in D_clusters.items():
        if len(v) != 0:
            D_clust_Dck[k] = {"Dck" : {}}
        else:
            pass

    # Get Dck values

    for ScafKey, cvalues in D_clusters.items():
        Get_Dckvalues(ScafKey, cvalues, DF_EvoDist, D_clust_Dck)

    ''' Get all the Dck values and calculate the DC statistic '''

    ALLscaffold_Dck = []

    # For each scaffold
    for Scaf, Values in D_clust_Dck.items():

        # Get the computed Dck values
        Scaf_Dck_values = list(Values["Dck"].values())

        # Add them to a common list
        ALLscaffold_Dck += Scaf_Dck_values

        # Compute the DC statistic, specific to a given scaffold
        # - Total num. of clusters across all scaffolds
        scf_allclusternum = len(Scaf_Dck_values)

        ''' DC Statistic - considering all the clusters in a given scaffold '''
        DCscf = DC_fx(scf_allclusternum, Scaf_Dck_values)
        D_clust_Dck[Scaf]["DC"] = DCscf

    ''' DC Statistic - considering all the clusters across the genome'''
    # Total number of clusters across all scaffolds
    ALLClustersNum = len(ALLscaffold_Dck)

    # DC glob
    DCglob = DC_fx(ALLClustersNum, ALLscaffold_Dck)

    # DT2 glob
    DT2glob = DT2_stat_glob(DF_EvoDist)

    # Cst glob
    CstGlob_stat = (DT2glob - DCglob)/DT2glob

    ''' DT1 Statistic - considering all the clusters in a given scaffold '''
    D_ScfClustFamGenes = {}

    for ScafKey, iClustInfo in D_clusters.items():
        if len(iClustInfo) != 0:
            D_ScfClustFamGenes[ScafKey] = {}

            # Variables to be filled
            ScafCgenenum = 0 # Num. of clustered genes in a scaffold
            ScafCgenelist = [] # List. of clustered genes in a scaffold

            for iclust in iClustInfo:
                # Get 
                # - cluster id
                # - number of clustered genes
                # - list of clustered genes
                Cname, Cgenenum, Cgenelist = iclust[0], iclust[2], iclust[3]

                # Discart overlapping genes
                discarted_genes = [i for i in Cgenelist if i.startswith("GG#") == True]

                if len(discarted_genes) != 0:
                    Cgenelist = [i for i in Cgenelist if i not in discarted_genes]
                    Cgenenum = Cgenenum - len(discarted_genes)
                    ScafCgenelist += Cgenelist

                else:
                    ScafCgenelist += Cgenelist

                # Save for later use
                for g in Cgenelist:
                    D_ScfClustFamGenes[ScafKey][g] = Cname

            # Save the list of clustered genes in this scaffold in a dictinoary (for later use with Mw test)
            D_ScfClustFamGenes[ScafKey]["AllClusteredGenes"] = ScafCgenelist

            # Number of clustered genes in this scaffold
            ScafCgenenum = len(ScafCgenelist)

            # Generate gene combinations
            ScafCgenecombo = list(combinations(ScafCgenelist, 2))

            # Get evol distances for each gene pair
            ScafCevodists = [DF_EvoDist.loc[j[0],j[1]] for j in ScafCgenecombo]

            # Compute DT1 for each scaffold
            DT1scf = DT_fx(ScafCgenenum, ScafCevodists)
            D_clust_Dck[ScafKey]["DT1"] = DT1scf
            # print("DTscf = ", DTscf)

        else:
            pass



    ''' DT2 Statistic - considering all the family genes in a given scaffold '''
    for ScafKey, iClustInfo in D_clusters.items():
        if len(iClustInfo) != 0:

            # Create a new key to access the dict with the gene list of family members in this scaffold
            tempKey = ScafKey.replace(".matrix","").split("_matrices.")[1]
            FamMembers = D_ScfFamGenes[tempKey]

            # Number of family genes in this scaffold
            ScafFCgenenum = len(FamMembers)

             # Generate gene combinations
            ScafFCgenecombo = list(combinations(FamMembers, 2))

            # Get evol distances for each gene pair
            ScafFCevodists = [DF_EvoDist.loc[j[0],j[1]] for j in ScafFCgenecombo]

            # Compute DT2 for each scaffold
            DT2scf = DT_fx(ScafFCgenenum, ScafFCevodists)
            D_clust_Dck[ScafKey]["DT2"] = DT2scf
            # print(ScafKey, "DT2scf = ", DT2scf)

    ''' Cst Statistic - considering all the family genes in a given scaffold '''

    for Scf, StatsInfo in D_clust_Dck.items():
        # print(Scf)
        DT2_stat = StatsInfo["DT2"]
        DC_stat = StatsInfo["DC"]

        Cst_stat = (DT2_stat - DC_stat)/DT2_stat
        StatsInfo["Cst"] = Cst_stat

        # for StatKey, StatValue in StatsInfo.items():
            # print("- ", StatKey, StatValue)


    ''' Prepare data for Mann-Whitney test'''

    MW_info = []

    for ScafKey, iClustInfo in D_clusters.items():
        if len(iClustInfo) != 0:

            # Create a new key to access the dict with the gene list of family members in this scaffold
            tempKey = ScafKey.replace(".matrix","").split("_matrices.")[1]
            FamMembers = D_ScfFamGenes[tempKey]

            # Number of family genes in this scaffold
            # ScafFCgenenum = len(FamMembers)

             # Generate gene combinations
            ScafFCgenecombo = list(combinations(FamMembers, 2))

            # Get evol distances for each gene pair
            # - Gene pairs are labeled as "Clustered" if both genes belong to the SAME cluster
            # - Gene pair "NotClustered" if both genes don't belong to any cluster

            for j in ScafFCgenecombo:
                gene1, gene2 = j[0], j[1]
                edist = DF_EvoDist.loc[gene1, gene2]

                if gene1 in D_ScfClustFamGenes[ScafKey]["AllClusteredGenes"] and gene2 in D_ScfClustFamGenes[ScafKey]["AllClusteredGenes"]:
                    g1_status = D_ScfClustFamGenes[ScafKey][gene1]
                    g2_status = D_ScfClustFamGenes[ScafKey][gene2]

                    if g1_status == g2_status:

                        info = [tempKey, gene1, gene2, g1_status, g2_status, "Clustered", edist]
                        MW_info.append(info)

                elif gene1 not in D_ScfClustFamGenes[ScafKey]["AllClusteredGenes"] and gene2 not in D_ScfClustFamGenes[ScafKey]["AllClusteredGenes"]:
                    info = [tempKey, gene1, gene2, "NC", "NC", "NotClustered", edist]
                    MW_info.append(info)
                else:
                    pass


    ''' Create a DataFrame with MW raw data '''
    DF_MWdata = pd.DataFrame.from_records(MW_info)
    DF_MWdata.columns = ["ScaffoldID", "Gene1ID", "Gene2ID", "Gene1_status", "Gene2_status", "GenePair_status", "EvoDist"]

    ''' Run Mann-Whitney test '''
    MW_test_results = []

    # Group by Scaffold
    DF_MWdata_gb = DF_MWdata.groupby("ScaffoldID")
    for s, sinfo in DF_MWdata_gb:

        # Empty dict
        temp_dict = {"Clustered" : 0, "NotClustered" : 0}

        # Fill it with count of gene pair status
        D_genepairstat = dict(Counter(sinfo["GenePair_status"]))
        temp_dict.update(D_genepairstat)

        # Get counts
        C_count = temp_dict["Clustered"]
        NC_count = temp_dict["NotClustered"]

        if C_count == 0 or NC_count == 0:
            # print("---",s, temp_dict)
            # MW test not possible
            pass
        else:
            # MW test possible
            # print("+++",s, temp_dict)

            NCvalues = sinfo[sinfo["GenePair_status"] == "NotClustered"]["EvoDist"]
            Cvalues = sinfo[sinfo["GenePair_status"] == "Clustered"]["EvoDist"]

            # Run MW test
            MWresults = Run_MannWhiteneyU(Cvalues, NCvalues)
            pvalue = MWresults[2]

            # Add Scaffold name
            MWresults.insert(0, s)

            # Add info about the sample size of "(Not)Clustered" points
            MWresults.insert(1, C_count)
            MWresults.insert(2, NC_count)

            # Add significance status
            MWresults.append(sigf_status(pvalue))
            MW_test_results.append(MWresults)

    # Create a dataframe
    DF_MWresults = pd.DataFrame.from_records(MW_test_results)
    DF_MWresults.columns = ["ScaffoldID", "ClusteredData", "NotClusteredData", "U1_Statistic", "U2_Statistic", "pvalue", "z_score", "MW_value", "variance_value", "stdev_value", "significance"]
    DF_MWresults["FamID"] = i_FAMname
    
    # Reorder Columns
    DF_MWresults = DF_MWresults.reindex(columns=["FamID", "ScaffoldID", "ClusteredData", "NotClusteredData", "U1_Statistic", "U2_Statistic", "pvalue", "significance", "z_score", "MW_value", "variance_value", "stdev_value"])

    # Add statistics
    for istat in ['DC', 'DT1', 'DT2', 'Cst']:
        add_stats(D_clust_Dck, DF_MWresults, istat)
    
    # Export results
    # Dck stats
    dck_outfile = f"{MWoutdir}/{i_FAMname}_Dck.rawdata.{i_gvhint}.csv"
    export_Dck(D_clust_Dck, dck_outfile)
    
    # MW raw data
    DF_MWdata["FamID"] = i_FAMname
    DF_MWdata = DF_MWdata.reindex(columns=["FamID", "ScaffoldID", "Gene1ID", "Gene2ID", "Gene1_status", "Gene2_status", "GenePair_status", "EvoDist"])

    DF_MWdata.to_csv(f"{MWoutdir}/{i_FAMname}_MannWhitney.rawdata.{i_gvhint}.csv", sep="\t", index=None)

    # MW results
    DF_MWresults.to_csv(f"{MWoutdir}/{i_FAMname}_MannWhitney.results.extended.{i_gvhint}.tsv", sep="\t", index=None)
    
    # MW results (the most important columns (for the report))
    DF_MWresults_brief = DF_MWresults.drop(columns=["U1_Statistic", "U2_Statistic", "z_score", "MW_value", "variance_value", "stdev_value", "DT1", "DT2", "DC"])
    DF_MWresults_brief = DF_MWresults_brief.reindex(columns=["FamID", "ScaffoldID","ClusteredData", "NotClusteredData", "Cst", "pvalue", "significance"])
    DF_MWresults_brief.columns = ["FamID", "ScaffoldID", "ClusteredGenes", "NotClusteredGenes", "Cst", "pvalue", "significance"]
    DF_MWresults_brief.to_csv(f"{MWoutdir}/{i_FAMname}_MannWhitney.results.brief.{i_gvhint}.tsv", sep="\t", index=None)
    
    return CstGlob_stat, DT2glob, DCglob, D_clust_Dck, DF_MWresults

# Get files
for kFAMILY, v in D_INPUTfiles.items():
    if "OverlappingGenesFile" not in v.keys():
        BEDfile, EvoDistfile, ClusterDict_list = v['BEDfile'], v['EvoDistfile'], v['ClusterDictfile']

        DF_FamBED = pd.read_csv(BEDfile, sep="\t", header=None)
        for ClusterDictfile in ClusterDict_list:
            gvalue_hint = get_gvalueID(ClusterDictfile)
    
            with open(f"{MWoutdir}/{kFAMILY}_GlobalStats_value.{gvalue_hint}.txt", "w") as o1:
                CstGlob_value, DT2glob_value, DCglob_value, b, c = GetStatistics_and_RunMWtest(kFAMILY, BEDfile, EvoDistfile, ClusterDictfile, gvalue_hint)
                print("FamID", "Statistic", "Value", sep="\t", end="\n", file=o1)
                print(kFAMILY, "Cst glob", CstGlob_value, sep="\t", end="\n", file=o1)
                print(kFAMILY, "Dt glob", DT2glob_value, sep="\t", end="\n", file=o1)
                print(kFAMILY, "Dc glob", DCglob_value, sep="\t", end="\n", file=o1)
    else:
        BEDfile, EvoDistfile, ClusterDict_list, OvlapFile = v['BEDfile'], v['EvoDistfile'], v['ClusterDictfile'], v['OverlappingGenesFile']

        DF_FamBED = pd.read_csv(BEDfile, sep="\t", header=None)
        for ClusterDictfile in ClusterDict_list:
            gvalue_hint = get_gvalueID(ClusterDictfile)

            with open(f"{MWoutdir}/{kFAMILY}_GlobalStats_value.{gvalue_hint}.txt", "w") as o1:
                CstGlob_value, DT2glob_value, DCglob_value, b, c = GetStatistics_and_RunMWtest(kFAMILY, BEDfile, EvoDistfile, ClusterDictfile, gvalue_hint, OvlapFile)        
                print("FamID", "Statistic", "Value", sep="\t", end="\n", file=o1)
                print(kFAMILY, "Cst glob", CstGlob_value, sep="\t", end="\n", file=o1)
                print(kFAMILY, "Dt glob", DT2glob_value, sep="\t", end="\n", file=o1)
                print(kFAMILY, "Dc glob", DCglob_value, sep="\t", end="\n", file=o1)


# In[ ]:

shutil.move(MWoutdir, clusterfinder_dir)



