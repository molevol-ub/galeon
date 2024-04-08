import ast, math, shutil, os, sys
from collections import Counter
import pandas as pd
from itertools import combinations
from scipy.stats import mannwhitneyu



# Input files
D_INPUTfiles = {}

# Search GFF(bed) files
I_GFF_list = [f"{sys.argv[1]}/{i}" for i in os.listdir(sys.argv[1]) if i.endswith("collapsed.temp")]

# Search distances matrices
I_distmatrices_list = [f"{sys.argv[2]}/{i}" for i in os.listdir(sys.argv[2]) if i.endswith("distancematrix.tsv")]

for r, d, l in os.walk(sys.argv[3]):
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
            
for r, d, l in os.walk(sys.argv[3]):
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

# Assign input files
for kFAM, v in D_INPUTfiles.items():
    for gfffile in I_GFF_list:
        if kFAM in gfffile:
            v["BEDfile"] = gfffile
            
    for evofile in I_distmatrices_list:
        if kFAM in evofile:
            v["EvoDistfile"] = evofile

rootname = os.path.abspath(sys.argv[3])
rootname = os.path.dirname(rootname)
MUoutdir = f"{rootname}/MannWhitney_StatisticsResults"

if os.path.exists(MUoutdir) == False:
    os.mkdir(MUoutdir)
else:
    shutil.rmtree(MUoutdir)
    os.mkdir(MUoutdir)


# In[4]:


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

    # Compute: 'mu', 'variance' and 'std. deviation'
    mu_value =  (n1*n2)/2 # (n1*n2)/2
    variance_value = (n1*n2*(n1+n2+1))/12
    stdev_value = math.sqrt(variance_value)

    # Compute 'z score'

    z_score = (abs(U_statistic - mu_value)-0.5)/stdev_value
    
    return z_score, mu_value, variance_value, stdev_value


def Run_MannWhiteneyU(x,y):
    U1_Statistic, pvalue = mannwhitneyu(x, y, alternative="two-sided")
    U2_Statistic = len(x)*len(y) - U1_Statistic

    z_score, mu_value, variance_value, stdev_value = other_statistics(len(x), len(y), U1_Statistic)
    
    return [U1_Statistic, U2_Statistic, pvalue, z_score, mu_value, variance_value, stdev_value]


# In[5]:


def GetStatistics_and_RunMUtest(i_FAMname, i_bedfile, i_evodistfile, i_clusterfile, i_gvalue):

    ''' Load GFF collapsed bed file '''
    # The genes are already sorted by the start coordinate

    DF_FamBED = pd.read_csv(i_bedfile, sep="\t", header=None)
    DF_FamBED.columns = ["Scaffold", "GeneStart", "GeneEnd", "Attribute", "bedclustID"]

    ''' Load Evolutionary Distance Matrix '''

    DF_EvoDist = pd.read_csv(i_evodistfile, sep="\t")
    DF_EvoDist = DF_EvoDist.set_index("Unnamed: 0")

    ''' Load Cluster dict '''

    CLUSTERdictfile = i_clusterfile

    with open(CLUSTERdictfile) as f1:
        temp = f1.read()

        D_clusters = ast.literal_eval(temp)


    # Empty dict to fill with the genelist IDs of each scaffold
    D_ScfFamGenes = {}

    # Group by scaffold and get the gene names

    DF_FamBED_gb = DF_FamBED.groupby("Scaffold")
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
                ScafCgenelist += Cgenelist

                # Save for later use
                for g in Cgenelist:
                    D_ScfClustFamGenes[ScafKey][g] = Cname

            # Save the list of clustered genes in this scaffold in a dictinoary (for later use with MU test)
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

    MU_info = []

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
                        MU_info.append(info)

                elif gene1 not in D_ScfClustFamGenes[ScafKey]["AllClusteredGenes"] and gene2 not in D_ScfClustFamGenes[ScafKey]["AllClusteredGenes"]:
                    info = [tempKey, gene1, gene2, "NC", "NC", "NotClustered", edist]
                    MU_info.append(info)
                else:
                    pass


    ''' Create a DataFrame with MU raw data '''
    DF_MUdata = pd.DataFrame.from_records(MU_info)
    DF_MUdata.columns = ["Scaffold", "Gene1", "Gene2", "Gene1_status", "Gene2_status", "GenePair_status", "EvoDist"]

    ''' Run Mann-Whitney test '''
    MU_test_results = []

    # Group by Scaffold
    DF_MUdata_gb = DF_MUdata.groupby("Scaffold")
    for s, sinfo in DF_MUdata_gb:

        # Empty dict
        temp_dict = {"Clustered" : 0, "NotClustered" : 0}

        # Fill it with count of gene pair status
        D_genepairstat = dict(Counter(sinfo["GenePair_status"]))
        temp_dict.update(D_genepairstat)

        # Get counts
        C_count = temp_dict["Clustered"]
        NC_count = temp_dict["NotClustered"]

        if C_count == 0 or NC_count == 0:
            pass
        else:
            # MU test possible
            NCvalues = sinfo[sinfo["GenePair_status"] == "NotClustered"]["EvoDist"]
            Cvalues = sinfo[sinfo["GenePair_status"] == "Clustered"]["EvoDist"]

            # Run MU test
            MUresults = Run_MannWhiteneyU(Cvalues, NCvalues)
            pvalue = MUresults[2]

            # Add Scaffold name
            MUresults.insert(0, s)

            # Add info about the sample size of "(Not)Clustered" points
            MUresults.insert(1, C_count)
            MUresults.insert(2, NC_count)

            # Add significance status
            MUresults.append(sigf_status(pvalue))
            MU_test_results.append(MUresults)

    # Create a dataframe
    DF_MUresults = pd.DataFrame.from_records(MU_test_results)
    DF_MUresults.columns = ["Scaffold", "Clustered Data", "NotClustered Data", "U1_Statistic", "U2_Statistic", "pvalue", "z_score", "mu_value", "variance_value", "stdev_value", "significance"]
    DF_MUresults["FamilyName"] = i_FAMname
    
    # Reorder Columns
    DF_MUresults = DF_MUresults.reindex(columns=["FamilyName", "Scaffold", "Clustered Data", "NotClustered Data", "U1_Statistic", "U2_Statistic", "pvalue", "significance", "z_score", "mu_value", "variance_value", "stdev_value"])

    # Export results
    # MU raw data
    DF_MUdata.to_csv(f"{MUoutdir}/{i_FAMname}_MannWhitney.rawdata_{i_gvalue}.csv", sep="\t", index=None)

    # MU results
    DF_MUresults.to_csv(f"{MUoutdir}/{i_FAMname}_MannWhitney.results_{i_gvalue}.tsv", sep="\t", index=None)


    with open(f"{MUoutdir}/{i_FAMname}_Stats_{i_gvalue}.dict", "w") as odict_file:
        odict_file.write(str(D_clust_Dck))
    


# Get files
for kFAMILY, v in D_INPUTfiles.items():
    BEDfile, EvoDistfile, ClusterDictfile = v['BEDfile'], v['EvoDistfile'], v['ClusterDictfile']
    DF_FamBED = pd.read_csv(BEDfile, sep="\t", header=None)

    for icfile in ClusterDictfile:
        temp = os.path.basename(icfile).split("temp_matrices_")[1]
        gvalue = temp.replace(".Dict.cluster.dict","")

        GetStatistics_and_RunMUtest(kFAMILY, BEDfile, EvoDistfile, icfile, gvalue)

