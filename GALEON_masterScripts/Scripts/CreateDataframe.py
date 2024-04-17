import pandas as pd
import os, ast, sys
from collections import Counter

''' Inputs '''
# Input coordinates file '''
CoordsFilePATH = sys.argv[1]

# Input cluster dictionary '''
DictFilePATH = sys.argv[2]

# Input chr size file '''
InputChrFile = sys.argv[3]

'''Output file '''
OutFile = sys.argv[4]

dirname = os.path.dirname(OutFile)

with open(DictFilePATH) as f1:
    temp2 = f1.read()
    D_Cluster = ast.literal_eval(temp2)
    

# Other parameters (Do not change them)
GFF3_attribute_field_sep = ";"
GFF3_geneID_start = "ID="
TwoFamAnalysis = False
D_TwoFamNames = {}

# Find the input coordinates file format
CoordsFileFormat = ""

for i in [".gff3", ".bed1", ".bed2"]:
    if i in CoordsFilePATH:
        CoordsFileFormat = i.replace(".","")
        
# checkpoint 1
if CoordsFileFormat == None:
    raise ValueError("Input Coordinate File has an unknown format")
else:
    pass


def get_geneIDs_for_Gff3(i_list):
    # Extract GeneIDs if the input coordinate file has a gff3 format
    Outinfo = []
    for i_attribute in i_list:
        temp = i_attribute.split(GFF3_attribute_field_sep)[0]
        geneID = temp.replace(GFF3_geneID_start,"")
        Outinfo.append(geneID)

    if len(Outinfo) == len(i_list):
        return Outinfo
    else:
        raise ValueError("Something is wrong, new attribute list length doesn'mat the original one")
        
def extract_features(i_scfname, i_list):
    # Extract features from Clustered genes information
    outrecords = []
    ClusterID, ClustereGenesNum, GeneNames = i_list[0], i_list[2], i_list[3]
    
    for idx, i in enumerate(GeneNames):
        outrecords.append([i_scfname, ClusterID, ClustereGenesNum, i, idx+1])
    return outrecords


def check_ovlap_genes(i_DF_ClustInfo):
    # Check for the presence of repeated genes
    temp = dict(Counter(i_DF_ClustInfo["GeneID"]))
    tempvalues = temp.values()
    repeated_genes = [k for k,v in temp.items() if v != 1]
    
    # get the list of overlapping clusters (corresponding to the repeated genes)
    ovlap_clustlist = []

    for k,v in i_DF_ClustInfo.groupby("GeneID"):
        if len(v["ClusterID"].unique()) != 1:
            temp = v["ClusterID"].unique()
            print(k, temp)

            for ii in temp:
                if ii not in ovlap_clustlist:
                    ovlap_clustlist.append(ii)
    
    # Export results
    if max(tempvalues) > 1:
        print(f"Warning! {len(repeated_genes)} Repeated geneIDs and {len(ovlap_clustlist)} Overlapping clusters were detected, this suggests the presence of overlapping clusters. The script will continue but we recommend to check the input data and consider trying a different g value")
        return True, repeated_genes, ovlap_clustlist # This will be considered in the Checkpoint 2
    else:
        return False, repeated_genes, ovlap_clustlist # This will be considered in the Checkpoint 2

''' Parse the input cluster dictionary '''
def ParseClusterDict(i_dict):
    # Save here all the records from the dictionary file
    OutRecords = []
    
    for ScfID, v in i_dict.items():
        if len(v) != 0:
            modKeyname = ScfID.split("_matrices.")[1].replace(".matrix","")
            for vi in v:
                info = extract_features(modKeyname, vi)
                OutRecords += info
    
    
    # Return a dataframe
    outDF = pd.DataFrame.from_records(OutRecords)
    try:
        outDF.columns = ["ScaffoldID", "ClusterID", "GenesNumber", "GeneID", "ClustID"]
        outDF = outDF.drop(["ClustID"], axis=1)
    except ValueError as e:
        outDF.columns = ["ScaffoldID", "ClusterID", "GenesNumber", "GeneID", "FamID", "ClustID"]
        outDF = outDF.drop(["ClustID"], axis=1)        
    
    return outDF

DF_ClustInfo = ParseClusterDict(D_Cluster)
OverlappingGenesStatus, rep_genes, ovlap_clusters = check_ovlap_genes(DF_ClustInfo)


# Load the input Coordinate file as dataframe
if CoordsFileFormat == "gff3" or CoordsFileFormat == "bed2" or CoordsFileFormat == "bed1":
    D_GeneCoords = {}
    
    # Load the GFF3
    DFcoords = pd.read_csv(CoordsFilePATH, sep="\t", header=None)
    
    # Number of input columns
    num_of_columns = len(DFcoords.columns)
    
    if num_of_columns == 5:
        # Name the columns
        DFcoords.columns = ["ScaffoldID", "GeneStart", "GeneEnd", "GeneID", "ClustID"]
        DFcoords = DFcoords.drop(["ClustID"], axis=1) # Drop ClusterID column which is useless

        # Replace the Attribute column
        DFcoords["GeneID"] = get_geneIDs_for_Gff3(DFcoords["GeneID"].to_list())
        
    elif num_of_columns == 6:
        # Name the columns
        DFcoords.columns = ["ScaffoldID",  "GeneStart", "GeneEnd", "GeneID", "FamID", "ClustID"]
        DFcoords = DFcoords.drop(["ClustID"], axis=1) # Drop ClusterID column which is useless
        # Replace the Attribute column
        DFcoords["GeneID"] = get_geneIDs_for_Gff3(DFcoords["GeneID"].to_list()) 
        
        TwoFamAnalysis = True # If True, then Two Families are being analysed at the same time and the script will work a bit different
        
    else:
        raise ValueError("Unknown error while loading coordinate file")
    
else:
    emsg = f"Error. Unknown 'CoordsFileFormat': {CoordsFileFormat}"
    raise ValueError(emsg)
# elif CoordsFileFormat in ["bed1", "bed2"]:
#     D_GeneCoords = {}
    
#     DFcoords = pd.read_csv(CoordsFilePATH, sep="\t", header=None)
#     DFcoords.columns = ["ScaffoldID", "GeneStart", "GeneEnd", "GeneID", "ClustID"]
#     DFcoords = DFcoords.drop(["ClustID"], axis=1)
# DFcoords



def CreateDF_for_UnclusteredGenes(i_DFClustInfo, i_DFCoords, iTwoFamAnalysis, oTwoFamdict):
    # Combine the information of Coordinates and Clustered genes to get the unclustered info
    outinfo = []
    
    if iTwoFamAnalysis == False:
        # Get Unclustered genes from the Coords file
        DF_unclustered = i_DFCoords[i_DFCoords["GeneID"].isin(i_DFClustInfo["GeneID"]) == False].copy()
        DF_unclustered_gb = DF_unclustered.groupby("ScaffoldID")

        for kScf, kDF in DF_unclustered_gb:
            for idx, (rowidx, row) in enumerate(kDF.iterrows()):
                info = [row["ScaffoldID"], f'#S{idx +1}', 1, row["GeneID"]]
                outinfo.append(info)
                
        # Create a DF
        outDF = pd.DataFrame.from_records(outinfo)

        # Name the columns
        outDF.columns = ["ScaffoldID", "ClusterID", "GenesNumber", "GeneID"]
        
        return outDF
    
    else:
        
        # If two families are processed at the same time, do this adjustment to parse the data
        i_DFCoords["GeneID_temp"] = i_DFCoords["FamID"] + "fam_" + i_DFCoords["GeneID"]
        
        ## Export those new gene names to dict
        
        for k,v in zip(i_DFCoords["GeneID"].to_list(), i_DFCoords["GeneID_temp"].to_list()):
            oTwoFamdict[k] = v
        
        # Get Unclustered genes from the Coords file
        DF_unclustered = i_DFCoords[i_DFCoords["GeneID_temp"].isin(i_DFClustInfo["GeneID"]) == False].copy()
        DF_unclustered_gb = DF_unclustered.groupby("ScaffoldID")

        for kScf, kDF in DF_unclustered_gb:
            for idx, (rowidx, row) in enumerate(kDF.iterrows()):
                info = [row["ScaffoldID"], f'#S{idx +1}', 1, row["GeneID"]]
                outinfo.append(info)
                
        # Create a DF
        outDF = pd.DataFrame.from_records(outinfo)

        # Name the columns
        outDF.columns = ["ScaffoldID", "ClusterID", "GenesNumber", "GeneID"]

        # Add an additional column with modified gene names (FamID + geneID)
        outDF["GeneID_temp"] = outDF["GeneID"].map(oTwoFamdict)    

        # Do the same with the other input dataframes
        D_TwoFamNames_r = {v:k for k,v in oTwoFamdict.items()}
        i_DFClustInfo["GeneID_temp"] = i_DFClustInfo["GeneID"].map(D_TwoFamNames_r)
        i_DFClustInfo.columns = list(i_DFClustInfo.columns[:-2]) + ["GeneID_temp", "GeneID"]

        return outDF

# get the information about the unclustered genes
DF_UnClust = CreateDF_for_UnclusteredGenes(DF_ClustInfo, DFcoords, TwoFamAnalysis, D_TwoFamNames)

# Checkpoint 2
if DF_ClustInfo.shape[0] + DF_UnClust.shape[0] != DFcoords.shape[0]:
    if OverlappingGenesStatus == True:
        DF_merged = pd.concat([DF_ClustInfo, DF_UnClust])
    else:
        emsg = f"Something is wrong. The sum of clustered ({DF_ClustInfo.shape[0]}) + unclustered ({DF_UnClust.shape[0]}) dataframe sizes, should match the total number of input genes from the Coordinate File: ({DFcoords.shape[0]})"
        raise ValueError(emsg)
else:
    DF_merged = pd.concat([DF_ClustInfo, DF_UnClust])


# DFcoords
if num_of_columns == 5:
    DF_Final = DF_merged.merge(DFcoords, how='inner', on='GeneID')
    DF_Final = DF_Final.drop(["ScaffoldID_y"], axis=1)
    DF_Final.columns = ["ScaffoldID", "ClusterID", "GenesNumber", "GeneID", "GeneStart", "GeneEnd"]
    
elif num_of_columns == 6:
    DF_Final = DF_merged.merge(DFcoords, how='inner', on='GeneID')
    DF_Final = DF_Final.drop(columns=["ScaffoldID_y", "GeneID_temp_y"], axis=1)
    DF_Final.columns = ['ScaffoldID', 'ClusterID', 'GenesNumber', 'GeneID', 'GeneID_original', 'GeneStart', 'GeneEnd', 'FamID']
else:
    raise ValueError("Unknown error")
# DF_Final


# Add chromosom sizes and names

with open(InputChrFile) as f1:
    D_ChrLength = {}
    D_ChrNames = {}
    
    if os.path.exists(InputChrFile):
        for line in f1:
            if line.startswith("#"):
                pass
            else:
                info = line.strip().split("\t")
                if len(info) == 3:
                    ScfID, ScfLength, ScfName = [*info]
                    D_ChrLength[ScfID] = ScfLength
                    D_ChrNames[ScfID] = ScfName
    else:
        raise ValueError("File not found!: Chromosome/Scaffold sizes")
        
    # Checkpoint 3
    if len(D_ChrLength) == 0 and len(D_ChrNames) == 0:
        raise ValueError("Error. Both dictionaries are empty!")
    elif len(D_ChrLength) != 0 and len(D_ChrNames) != 0:
        if len(D_ChrLength) != len(D_ChrNames):
            raise ValueError("Error. Dictionaries' size doesn't match!")
        else:
            # print(len(D_ChrLength), len(D_ChrNames))
            
            # Add scaffold lengths
            DF_Final["ScfLength"] = DF_Final["ScaffoldID"].map(D_ChrLength)
            
            # Add scaffd names
            DF_Final["ScfName"] = DF_Final["ScaffoldID"].map(D_ChrNames)
    
    elif len(D_ChrLength) != 0 and len(D_ChrNames) == 0:
        DF_Final["ScfLength"] = DF_Final["ScaffoldID"].map(D_ChrLength)
    
    else:
        pass
        
# Change the order of the output columns
if num_of_columns == 5:
    newColOrder = ["ScaffoldID", "ScfLength", "ScfName", "ClusterID", "GenesNumber", "GeneID", "GeneStart", "GeneEnd"]
    DF_Final = DF_Final[newColOrder]
    
elif num_of_columns == 6:
    newColOrder = ['ScaffoldID', 'ScfLength', 'ScfName',
                   'ClusterID', 'GenesNumber', 'FamID','GeneID_original', 'GeneID',
                   'GeneStart', 'GeneEnd']
    DF_Final = DF_Final[newColOrder]
else:
    raise ValueError("Unknown error at the final step of the script")


def Twofam_composition(i_famlist, i_genelist):
    # Open keys
    D_counter = {i : 0 for i in i_famlist}

    for geneID in i_genelist:
        famID = geneID.split("fam_")[0]
        D_counter[famID] += 1

    tlist = list(D_counter.values())
    twofam_clustinfo = f"{sum(tlist)} ({tlist[0]}, {tlist[1]})"
    return twofam_clustinfo

def AddTables(i_DF, o_filename, i_opt=False, i_emptyDF=False):

    # drop some columns
    if i_opt == True:
        tempDF = i_DF.drop(columns = ["GeneID", "GeneStart", "GeneEnd", "GeneID_original"])
    else:
        tempDF = i_DF.drop(columns = ["GeneID", "GeneStart", "GeneEnd"])

    # drop duplicates
    dropdups = tempDF.drop_duplicates()

    # Sort by scaffold length
    dropdups = dropdups.sort_values(by="ScfLength", ascending=False)

    # Export the table with cluster sizes
    o_file = o_filename.replace("GeneLocation", "ClusterSizes")
    dropdups.to_csv(o_file, index=None, sep="\t")

    ###

    # Create a new dataframe with a summary of clustered/singleton genes
    temp_list = []
    for i in dropdups["ClusterID"].to_list():
        if i.startswith("#S"):
            temp_list.append("Singleton")
        else:
            temp_list.append("Clustered")
    dropdups["Category"] = temp_list

    outlist = []
    for k,v in dropdups.groupby(["ScaffoldID", "Category"]):
        scfID, chrID = k[0], k[1]
        if i_opt == True:
            info = [scfID, v["ScfLength"].to_list()[0], v["ScfName"].to_list()[0], chrID, sum(v["GenesNumber"].to_list()), v["FamID"].to_list()[0]]
        else:
            info = [scfID, v["ScfLength"].to_list()[0], v["ScfName"].to_list()[0], chrID, sum(v["GenesNumber"].to_list()), v["FamID"].to_list()[0]]
        outlist.append(info)

    # Summary DF
    SumDF = pd.DataFrame.from_records(outlist)

    SumDF.columns = ["ScaffoldID", "ScfLength", "ScfName", "Category", "GenesNumber", "FamID"]
    SumDF = SumDF.reindex(columns=["FamID", "ScaffoldID", "ScfLength", "ScfName", "Category", "GenesNumber"])

    # if i_opt == True:
    #     SumDF.columns = ["ScaffoldID", "ScfLength", "ScfName", "Category", "GenesNumber", "FamID"]
    #     SumDF = SumDF.reindex(columns=["FamID", "ScaffoldID", "ScfLength", "ScfName", "Category", "GenesNumber"])
    # else:
    #     SumDF.columns = ["ScaffoldID", "ScfLength", "ScfName", "Category", "GenesNumber", "FamID"]
    #     SumDF = SumDF.reindex(columns=["FamID", "ScaffoldID", "ScfLength", "ScfName", "Category", "GenesNumber"])

    SumDF = SumDF.sort_values(by=["ScfLength","Category"], ascending=False)

    # Export the table with gene organization info
    o_file = o_filename.replace("GeneLocation", "GeneOrganizationSummary")
    SumDF.to_csv(o_file, index=None, sep="\t")

    # Export the table with cluster sizes over the entire genome
    o_file = o_filename.replace("GeneLocation", "GeneOrganizationGenomeSummary")

    tempDF = SumDF[["FamID","Category","GenesNumber"]].copy(deep=True)
    temprecs = []
    for k,v in tempDF.groupby("Category"):
        if OverlappingGenesStatus == False:
            info = [v["FamID"].unique()[0], k, sum(v["GenesNumber"])]
            temprecs.append(info)

        elif OverlappingGenesStatus == True:
            if k == "Clustered":
                info = [v["FamID"].unique()[0], k, sum(v["GenesNumber"])]
                temprecs.append(info)

                info = [v["FamID"].unique()[0], "Repeated genes", len(rep_genes)]
                temprecs.append(info)

                info = [v["FamID"].unique()[0], "Overlapping clusters", len(ovlap_clusters)]
                temprecs.append(info)

            else:
                info = [v["FamID"].unique()[0], k, sum(v["GenesNumber"])]
                temprecs.append(info)

    SumDF_glob = pd.DataFrame.from_records(temprecs)
    SumDF_glob.columns = ["FamID","Category","GenesNumber"]
    SumDF_glob.to_csv(o_file, index=None, sep="\t")

# Export raw data results
DF_Final["ScfLength"] = DF_Final["ScfLength"].apply(pd.to_numeric)
if dirname.startswith("merged_"):
    pass
else:
    DF_Final["FamID"] = dirname

# reindex
if TwoFamAnalysis == True:
    newcolorder = ['FamID', 'ScaffoldID', 'ScfLength', 'ScfName', 'ClusterID', 'GenesNumber', 'GeneID_original', 'GeneID', 'GeneStart', 'GeneEnd']
elif TwoFamAnalysis == False:
    newcolorder = ['FamID', 'ScaffoldID', 'ScfLength', 'ScfName', 'ClusterID', 'GenesNumber', 'GeneID', 'GeneStart', 'GeneEnd']

DF_Final = DF_Final.reindex(columns=newcolorder)

DF_Final = DF_Final.sort_values(by="ScfLength", ascending=False)
DF_Final.to_csv(OutFile, sep="\t", index=None)

# Add more tables
if TwoFamAnalysis == False:
    AddTables(DF_Final, OutFile, TwoFamAnalysis)
else:
    ####
    # Group by ScfID and Cluster ID
    # Retain only the info of clusters with members of both families
    outlist = []
    twofam_clusterskeys = []
    twofamID = ""

    for k,v in DF_Final.groupby(["ScaffoldID", "ClusterID"]):
        if len(v["FamID"].unique()) != 1:
            famnames_sorted = sorted(v["FamID"].unique())
            iScfID = k[0]
            iScfLen = v["ScfLength"].unique()[0]
            iChrID = v["ScfName"].unique()[0]
            iClusterID = k[1]
            iClustSize = v["GenesNumber"].unique()[0]
            iClusterComposition = Twofam_composition(v["FamID"].unique(), v["GeneID"].to_list())
            iFamName = " and ".join(famnames_sorted)
            twofamID = ".".join(v["FamID"].unique())
            
            info = [iScfID, iScfLen, iChrID, iClusterID, iClustSize, iClusterComposition, iFamName]
            outlist.append(info)
            
            twofam_clusterskeys.append([iScfID, iClusterID])
            
    # create the summary DF
    if len(twofam_clusterskeys) != 0:
        mDF = pd.DataFrame.from_records(outlist)
        mDF.columns = ["ScaffoldID", "ScfLength", "ScfName", "ClusterID", "GenesNumber", "Composition", "FamID"]

        # reindex
        mDF = mDF.reindex(columns=["FamID", "ScaffoldID", "ScfLength", "ScfName", "ClusterID", "GenesNumber", "Composition"])
        
        # Out name
        filename = os.path.basename(OutFile).replace("merged_family",f"{twofamID}_family").replace("GeneLocation", "GeneOrganizationSummary")
        OutFile2 = f"{dirname}/{filename}"
        mDF.to_csv(OutFile2, sep="\t", index=None)
    else:
        pass

    ### --- end ---

    ####
    # Create a new dataframe by excluding all the clusters where two fam. members are present
    if len(twofam_clusterskeys) != 0:
        filtered_list = []
        for idx, row in DF_Final.iterrows():
            info = row.loc[["ScaffoldID","ClusterID"]].to_list()
            if info in twofam_clusterskeys:
                pass
            else:
                filtered_list.append(row)

        DF_oneFam = pd.DataFrame.from_records(filtered_list)
    else:
        DF_oneFam = DF_Final

    # checkpoint 4: now only clusters formed by one family shoud be present
    for k,v in DF_oneFam.groupby(["ScaffoldID", "ClusterID"]):
        famNum = len(v["FamID"].unique())# num of families present in a give cluster
        if famNum != 1:
            emsg = f"Two families present in the same cluster: {k}, something is wrong here"
            raise ValueError(emsg)
            
    #### --- end ---

    ####
    famlist =  list(DF_Final["FamID"].unique())
    for ifam in famlist:
        tempDF = DF_oneFam[DF_oneFam["FamID"] == ifam].copy(deep=True)

        filename = os.path.basename(OutFile).replace("merged_family",f"{ifam}_family")
        OutFile3 = f"{dirname}/{filename}"

        if len(tempDF) != 0:
            AddTables(tempDF, OutFile3, TwoFamAnalysis)

        else:
            # Export empty tables
            
            columnlist = tempDF.columns
            x = ["-"]*(len(columnlist)-1)
            x.insert(0,[ifam])

            eDF = pd.DataFrame.from_records(x).T
            eDF.columns = columnlist

            ####
            tempDF_E1 = eDF.drop(columns = ["GeneID", "GeneStart", "GeneEnd", "GeneID_original"])
            o_file = OutFile3.replace("GeneLocation", "ClusterSizes")
            tempDF_E1.to_csv(o_file, index=None, sep="\t")

            ####
            x = ["-"]*5
            x.insert(0,[ifam])
            eDF2 = pd.DataFrame.from_records(x).T
            eDF2.columns = ["FamID", "ScaffoldID", "ScfLength", "ScfName", "Category", "GenesNumber"]
            o_file = OutFile3.replace("GeneLocation", "GeneOrganizationSummary")
            eDF2.to_csv(o_file, index=None, sep="\t")

            ####
            x = ["-"]*2
            x.insert(0,[ifam])
            eDF3 = pd.DataFrame.from_records(x).T
            eDF3.columns = ["FamID","Category","GenesNumber"]
            o_file = OutFile3.replace("GeneLocation", "GeneOrganizationGenomeSummary")
            eDF3.to_csv(o_file, index=None, sep="\t")            
    ### --- end ---
