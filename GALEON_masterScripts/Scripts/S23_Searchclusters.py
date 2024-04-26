import sys, os, time, shutil
import pandas as pd
import copy as cp
from itertools import combinations

''' Parameteres '''


g_param_value = float(sys.argv[1]) # "100"

ExpandControl_param = sys.argv[2] # inf o integer (1,2,3) "1"

if ExpandControl_param == "inf":
    pass
else:
    try:
        ExpandControl_param = int(ExpandControl_param)
    except ValueError:
        raise ValueError("Error in 'ExpandControl_param'. Acceptable values are: 'inf' or an integer value (1,2,3,etc) ")

directory = sys.argv[3]

FamNumber_param = sys.argv[4] # WithinFamilies / BetweenFamilies

''' Functions '''
## Parametres de la matriu
def i_mx_param(i_matrix):
    genenum = i_matrix.shape[0]
    return list(range(0,genenum))

## Com es definiran els clústers?
def CL_criteria(g_parameter, num_of_genes): 
    # g_parameter => 1 unitat Cl
    # num_of_genes => maximum number of genes that may form a cluster (the largest theoretical cluster will include all the genes from a given scaffold)
    # els clústers més grans. Es pot posar qualsevol valor
    
    Cl_dict = {}
    gnum_range = list(range(num_of_genes+1))

    for n_i in gnum_range[2:]:
        Cl_value = g_parameter*(n_i-1)
        Cl_dict[n_i] = Cl_value
    return Cl_dict


# Detectar els potencials clústers
def Find_Clusters(i_mx, i_mx_range, i_mx_max, cl_dictionary):

    defined_clusters  = {}

    for i in i_mx_range:
        n = 1
        start = i
        fw_step = i+1
        # end = [start]
        print(f"## Start iteration from '{start}' index")
        defined_clusters[i] = [start]

        while fw_step <= i_mx_max:
            # Distàncies inici, pas
            n_plus = n + 1
            dist_one = i_mx[i,i]

            dist_fw_step = i_mx[start,fw_step]
            fw_step += 1
            # print(dist_fw_step)
            cl_threshold = cl_dictionary[n_plus]
            # print("!",cl_threshold)

            if dist_fw_step <= cl_threshold:
                # print(True)
                n = n_plus
                # end.append(fw_step-1)
                defined_clusters[i].append(fw_step-1)
                print(f"{start} : {start}->{fw_step-1} = {dist_fw_step} => Yes, n = {n}, Cl = {cl_threshold}")
            
            else:
                n = n_plus

                print(f"{start} : {start}->{fw_step-1} = {dist_fw_step} => No, n = {n}, Cl = {cl_threshold}")
                print(f"## -- End of cluster search from '{start}' -- ##\n" )
                break

        else:
            print(f"## -- End of cluster search from '{start}' -- ##\n" )
            pass

    defined_clusters_raw = {k:v for (k,v) in defined_clusters.items() if len(v) > 1}
    print("## Matrix processed ##")
    print("## Job done ##\n")
    return defined_clusters_raw


## Descartar clusters niats (uns dins dels altres) => intersecció 100%
## Ex: cluster 1 [0,1,2,3], cluster 2 [1,2,3] => el cluster 1 conté el 2 i per tant descarto el 2

class Check_data:

    def __init__(self, w_tuple):
        self.w_tuple = w_tuple

    def return_dict_values(self, i_tuple, i_dict): # dict original
        cset1, cset2 = i_dict[i_tuple[0]], i_dict[i_tuple[1]]

        return cset1, cset2
    
def remove_dict_keys(i_klist, i_dict_dcp): # deepcopy
    for ikr in i_klist:
        del i_dict_dcp[ikr]

def cluster_redundancy(clust_rawdict):
    clust_rawdict_nr = cp.deepcopy(clust_rawdict)
    idx_comb = list(combinations(clust_rawdict.keys(),2)) # combinacio de claus del diccionari
    if len(idx_comb) > 0:
        # idx_comb_filt = [ii for ii in idx_comb if ii[1] == ii[0] + 1] # em quedo amb parelles de clusters en tàndem
        # per tal de poder veure si hi ha clusters solapats o no
        k_to_remove = set()
        print("-- Matrix index combinations",idx_comb)
        for cidx in idx_comb:
            cidx_inst = Check_data(cidx)
            # c_set1 = set()
            # c_set2 = set()

            # if cidx_inst.check_in_list(k_to_remove) == False: # si no hi ha intersecció entre els descartats i els evaluats
            csets = cidx_inst.return_dict_values(cidx,clust_rawdict)
            c_set1, c_set2 = set(csets[0]), set(csets[1])

            if c_set1.issuperset(c_set2): # Si surt == True => c_set2 està completament dins de c_set1
                    # k_to_remove.append(cidx[1])
                    # del clust_rawdict[cidx[1]]
                k_to_remove.add(cidx[1])
                print(f"! Complete overlap found between {cidx}")#=> removing '{cidx[1]}' info from input dictionary")


            elif c_set2.issuperset(c_set1): # Si surt == True => c_set1 està completament dins de c_set2
                    # k_to_remove.append(cidx[0])
                    # del clust_rawdict[cidx[0]]
                k_to_remove.add(cidx[0])
                print(f"! Complete overlap found between {cidx}")#" => removing '{cidx[0]}' info from input dictionary")        

            else:
                pass

        print("## - removing redundant information - ##")

        print(k_to_remove)

        remove_dict_keys(k_to_remove,clust_rawdict_nr)

        print("## - done - ##")

        return clust_rawdict_nr
    else:
        print("## - No redundant clusters have been found as there is only one cluster in this scaffold - ##")
        print("## - done - ##")
        return(clust_rawdict)

# print("## (2) ## Checkpoint: looking for nested clusters ###########")

# mx_clusters_nr = cluster_redundancy(mx_clusters)
# print(mx_clusters_nr)
# print("\n")

# print("## (3) ## Checkpoint: looking for overlapping clusters  ###########")

def cluster_overlap(clust_dict_nr):
    idx_comb = list(combinations(clust_dict_nr.keys(),2)) # combinacio de claus del diccionari
    # per tal de poder veure si hi ha clusters solapats o no
    print("-- Matrix index combinations",idx_comb)

    overlap_list = []
    for cidx in idx_comb:
        c_set1, c_set2 = set(clust_dict_nr[cidx[0]]), set(clust_dict_nr[cidx[1]])
        # overlap_gindexes = []
        # if c_set1.isdisjoint(c_set2) == False:
        if c_set1.isdisjoint(c_set2) == False:
            overlap_gindexes = c_set1.intersection(c_set2)
            print(f"! Partial overlapping found between {cidx}")
            cidx_pair = f"{cidx[0]}_{cidx[1]}"
            overlap_list.append([cidx_pair,list(overlap_gindexes)])

            # print(c_olap_geneids)
        else:
            pass
            
    olap_list_len = len(overlap_list)
    
    if olap_list_len != 0:
        print(f"-- 'Overlapping clusters' cases: {olap_list_len} --")

    return overlap_list

def rename_expanded_cl_dict(ii_dict):
    newkeys = [min(vi) for vi in ii_dict.values()]
    r_ii_dict = dict(zip(newkeys,ii_dict.values()))
    return r_ii_dict

def cluster_expansion(ii_mx, clust_dict_1round, ii_cl_dictionary, i_totgeneinscf):
    print("====(1)====>",clust_dict_1round)

    # idx_combi = list(combinations(clust_dict_1round.keys(),2))
    # print("@~@~@~@~@~@~@~@~@", clust_dict_1round)

    # print(clust_dict_1round)
    for k,v in clust_dict_1round.items():
        start = k
        pot_clust_genenum = len(v) + 1 # el nou número de potencial de gens en clústers = gens_encluster + 1 (va incrementant de 1 en 1 si es pot ampliar el clúster)
        
        oneback_step = k-1 # partint des del primer gen del clúster, es torna un gen enrere per veure si es pot ampliar
        if oneback_step >= 0:
            last_cl_gene = v[-1] # des del gen anterior, es mira la distància fins l'últim gen del clúster d'entrada. 
            # Ex: cluster entre gens [3,4] => el primer gen és el 3. Llavors, tornaríem al gen 2 i miraríem la distància fins el 4 <=> de 2 -> 4

            dist_back_step = ii_mx[oneback_step, last_cl_gene]

            try:
                cl_threshold_exp = ii_cl_dictionary[pot_clust_genenum]
                # print(f"Start: {k}")
                # print(f"Backstep: {oneback_step}")
                print(f"- Distance to the last cluster gene ({oneback_step}->{last_cl_gene}): {dist_back_step}")
                
                while dist_back_step <= cl_threshold_exp and oneback_step >= 0: # si la dist cau per sota o és igual al llindar establert per Cl, s'amplia el clúster
                    print(f"\tTrue: {dist_back_step} <= {cl_threshold_exp}, +1 gene to {start}_id cluster")


                #     # pot_clust_genenum += 1 # si es compleix la condició, s'amplia ampliar el clúster => e
                    # print(v)

                    clust_dict_1round[k].insert(0,oneback_step) # si es compleix la condició, s'amplia ampliar el clúster => e
                    pot_clust_genenum += 1
                    # print(v)

                    oneback_step -= 1 # tornem un gen enrere (si s'ha pogut ampliar el clúster abans)
                    if oneback_step < 0:
                        pass
                    else:
                        # last_cl_gene -= 1

                        dist_back_step = ii_mx[oneback_step, last_cl_gene]
                        cl_threshold_exp = ii_cl_dictionary[pot_clust_genenum]

                        print(f"\t-- new expansion test in progress ({oneback_step}->{last_cl_gene}) --")

                        if dist_back_step <= cl_threshold_exp:
                            print("\t--> one more added --\n")
                        else:
                            print(f"\n\tFalse: {dist_back_step} > {cl_threshold_exp}, no gene added to {start}_id cluster")
                            print("\t-- end --")
                else:
                    pass
            except KeyError:
                if pot_clust_genenum > i_totgeneinscf:
                    print("# Reached the first gene of the scaffold. Nothing else can be done here.")
                    pass
        else:
            print("# No expansion as the beggining of the matrix was reached - ##")
            pass

    clust_dict_2round = rename_expanded_cl_dict(clust_dict_1round)
    print("====(2)====>",clust_dict_2round)

    return clust_dict_2round


def export_results(clust_dict_ref, geneindex_dict, ofile, opt, o_dict, df_filename, i_ovlap_list):

    klist = list(clust_dict_ref.keys())
    vlist = list(clust_dict_ref.values())
    vgene = []

    for vii in vlist:
        # try:
        vii_list = [geneindex_dict[vid] for vid in vii]
        vgene.append(vii_list)
        # except KeyError:
            # break

    kidx_dict = { k1 : [f"Cluster{klist.index(k1) + 1}", f"matrix_ref_key={k1}", len(clust_dict_ref[k1]), vgene[klist.index(k1)]] for k1 in klist}
    
    for k2,v in clust_dict_ref.items():
        clustinfo, genesinclust = kidx_dict[k2][0], kidx_dict[k2][2]
        cinfo = f">{clustinfo}\t{genesinclust}\n"
        if opt == 1:
            ofile.write(cinfo)
            for i in v:
                ofile.write(f"--{geneindex_dict[i]}\n") 
        elif opt ==    2:
            ofile.write(cinfo)
            oinfo = kidx_dict[k2]
            if oinfo not in o_dict[df_filename]:
                o_dict[df_filename].append(oinfo)
                print(oinfo)


    if opt == 3:
        ofile.write("#################\n")
        ofile.write(df_filename+"\n")
        for oi in i_ovlap_list:
            case_num = i_ovlap_list.index(oi) + 1
            oi_keys = oi[0].split("_")
            oi_values = [geneindex_dict[iii] for iii in oi[1]]
            ofile.write(f"Overlap case {case_num}:\n")
            ofile.write(f"{kidx_dict[int(oi_keys[0])][0]} and {kidx_dict[int(oi_keys[1])][0]}: {oi_values}\n\n")


def process_matrices(i_matrix_list, o_explicitfile, o_shortfile, o_overlapfile, o_dictfile, location, g_param, exp_control_param):
    outdict = {}
    for df_file in i_matrix_list:
        outdict[df_file] = []
        filen = f"\n# INPUT: {df_file}\n" 
        print(filen)

        with open(f"{location}/{df_file}") as f:
            # import_file = pd.read
            pdmatrix = pd.read_csv(f,sep="\t")
            pd.set_option("display.max_rows",5)
            pd.set_option("display.max_columns",5)

            gene_index_assoc = {}
            genelist = list(pdmatrix.columns[1:])
            genenumber = len(genelist) # total genes in a given scaffold
            col_name = pdmatrix.columns[0]
            geneindex = range(0,genenumber)

            gene_index_assoc = {gindex : gname for (gindex, gname) in zip(geneindex,genelist)}
            pdmatrix.drop(col_name,axis=1,inplace=True)
            npmatrix = pdmatrix.to_numpy()
            npmatrix_t = npmatrix.transpose()

            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
            print(f"{g_param}kb, {g_param*1e3}bp")
            # npmatrix_t = npmatrix_t/(g_param*1e3) # Ex: input 100 g value = 100 kb => Matrix must be divided by 100*1e3 = 100_000 value (bp units) # update (22 Abril 2024)
            npmatrix_t = npmatrix_t/1e3 # Ex: input 100 g value = 100 kb => Matrix (in bp) must be divided by 1e3 => now is scaled to kb # update (26 Abril 2024)


            print("## (0) ## Matrix pre-processing  ###########")

            mx_range = i_mx_param(npmatrix_t) # rang de valors (que serviran per a moure's per moure's sobre la matriu)
            mx_max = max(mx_range)


            print(f"## - Clustering criteria: CL = {g_param} ; Gene number: {genenumber} by default ## ")
            mx_Cl = CL_criteria(g_param, genenumber) ### arguments ### default: 1 (1g), total number of genes in a scaffold (thus, the largest scaffold should include all the genes at most)

            print("## (1) ## Detection of potential Clusters ###########")

            mx_clusters = Find_Clusters(npmatrix_t, mx_range, mx_max, mx_Cl)

            if len(mx_clusters) != 0:
            
                o_explicitfile.write(filen)
                o_shortfile.write(filen)
            
                print("## (1.1) ## Checkpoint: looking for nested clusters ###########")

                mx_clusters_nr = cluster_redundancy(mx_clusters)
                print("\n")

                print("## (1.2) ## Checkpoint: looking for overlapping clusters  ###########")

                first_overlap_test = cluster_overlap(mx_clusters_nr)
                print("\n")

                # As far as it is known, one round of gene cluster search and one expansion round are enough.
                # However, just in case more than one expansion round is possible, a while loop has been implemented.
                # If exp_control_param is <1, that is, if it is '0', no Expansion round is done.

                Expand_control = True 
                Expand_control_count = 1

                if exp_control_param == "inf":
                    print("## (2) ## Trying to Expand 1st Round Clusters - ENABLED ###########")
                
                elif Expand_control_count <= exp_control_param:
                    print("## (2) ## Trying to Expand 1st Round Clusters - ENABLED ###########")
                else:
                    print("## (2) ## Trying to Expand 1st Round Clusters - DISABLED ! ###########")

                # Unlimited expansion
                if exp_control_param == "inf":
                    while Expand_control == True:

                        print(f"EXPAND ITERATION: {Expand_control_count}")
                        clusters_def_2round = cluster_expansion(npmatrix_t, mx_clusters_nr, mx_Cl, genenumber)
                        print("\n")

                        if mx_clusters_nr == clusters_def_2round:

                            print("## ! (2) ## Checkpoint: looking for nested clusters ###########")

                            clusters_def_2round_nr = cluster_redundancy(clusters_def_2round)

                            print("\n")

                            print("## ! (2) ## Checkpoint: looking for overlapping clusters  ###########")
                            last_overlap_test = cluster_overlap(clusters_def_2round_nr)
                            print("\n")

                            print(f"- Stop expanding after '{Expand_control_count}' iteration ")
                            Expand_control = False

                            print("## (3) ## Exporting results  ###########")

                            export_results(clusters_def_2round_nr,gene_index_assoc,o_explicitfile,1,"","","")
                            export_results(clusters_def_2round_nr,gene_index_assoc,o_shortfile,2,outdict,df_file,"")

                            if len(last_overlap_test) != 0: 
                                export_results(clusters_def_2round_nr,gene_index_assoc,o_overlapfile,3,outdict,df_file,last_overlap_test)

                            else:
                                pass

                        else:
                            Expand_control_count += 1 # update Expand_control_count variable value

                            print("## ! (2) ## Checkpoint: looking for nested clusters ###########")

                            clusters_def_2round_nr = cluster_redundancy(clusters_def_2round)

                            print("\n")

                            print("## ! (2) ## Checkpoint: looking for overlapping clusters  ###########")
                            last_overlap_test = cluster_overlap(clusters_def_2round_nr)
                            print("\n")

                            print("## (3) ## Exporting results  ###########")

                            export_results(clusters_def_2round_nr,gene_index_assoc,o_explicitfile,1,"","","")
                            export_results(clusters_def_2round_nr,gene_index_assoc,o_shortfile,2,outdict,df_file,"")

                            mx_clusters_nr = clusters_def_2round_nr

                            if len(last_overlap_test) != 0: #???????????????????? ficar-ho en el loop
                                export_results(clusters_def_2round_nr,gene_index_assoc,o_overlapfile,3,outdict,df_file,last_overlap_test)

                            else:
                                pass
                
                # Expansion limited to user defined parameter 
                
                else:
                    while Expand_control == True and Expand_control_count <= exp_control_param:

                        print(f"EXPAND ITERATION: {Expand_control_count} out of {exp_control_param}")
                        clusters_def_2round = cluster_expansion(npmatrix_t, mx_clusters_nr, mx_Cl, genenumber)
                        print("\n")

                        if mx_clusters_nr == clusters_def_2round:

                            print("## ! (2) ## Checkpoint: looking for nested clusters ###########")

                            clusters_def_2round_nr = cluster_redundancy(clusters_def_2round)

                            print("\n")

                            print("## ! (2) ## Checkpoint: looking for overlapping clusters  ###########")
                            last_overlap_test = cluster_overlap(clusters_def_2round_nr)
                            print("\n")

                            print(f"- Stop expanding after '{Expand_control_count}' iteration out of {exp_control_param}")
                            Expand_control = False

                            print("## (3) ## Exporting results  ###########")

                            export_results(clusters_def_2round_nr,gene_index_assoc,o_explicitfile,1,"","","")
                            export_results(clusters_def_2round_nr,gene_index_assoc,o_shortfile,2,outdict,df_file,"")

                            if len(last_overlap_test) != 0: #???????????????????? ficar-ho en el loop
                                export_results(clusters_def_2round_nr,gene_index_assoc,o_overlapfile,3,outdict,df_file,last_overlap_test)

                            else:
                                pass

                        else:
                            Expand_control_count += 1 # update Expand_control_count variable value

                            print("## ! (2) ## Checkpoint: looking for nested clusters ###########")

                            clusters_def_2round_nr = cluster_redundancy(clusters_def_2round)

                            print("\n")

                            print("## ! (2) ## Checkpoint: looking for overlapping clusters  ###########")
                            last_overlap_test = cluster_overlap(clusters_def_2round_nr)
                            print("\n")

                            print("## (3) ## Exporting results  ###########")

                            export_results(clusters_def_2round_nr,gene_index_assoc,o_explicitfile,1,"","","")
                            export_results(clusters_def_2round_nr,gene_index_assoc,o_shortfile,2,outdict,df_file,"")

                            mx_clusters_nr = clusters_def_2round_nr

                            if len(last_overlap_test) != 0: #???????????????????? ficar-ho en el loop
                                export_results(clusters_def_2round_nr,gene_index_assoc,o_overlapfile,3,outdict,df_file,last_overlap_test)

                            else:
                                pass
                    
                
            else:
                print(f"# No clusted were identified in this scaffold, for this 'g value' = {g_param_value}")

            print("\n## Job done :) ##\n\n")

    o_dictfile.write(str(outdict))

''' Execution '''
    
start = time.perf_counter()

outdict = {}

wkd_dirlist = [i[0] for i in os.walk(directory)][1:]
wkd_dirlist = [i for i in wkd_dirlist if i.endswith("_matrices")]

if len(wkd_dirlist) == 0:
    msg = f"No input directories (expected directory format: *_matrices) have been found in this wkd: {directory} !!!"
    raise ValueError(msg)

for subfolder in wkd_dirlist:
    print("## INPUT directory: ",subfolder)

    # Create a specific directory for each the test g_value 
    # gparam_shortID = str(g_param_value / 100) # g param short ID
    gparam_shortID = str(g_param_value) # g param short ID (update, 22 Abril 2024)
    subfoldername = "/".join(subfolder.split("/")[0:])

    gparam_specific_dir = f"{subfoldername}_{gparam_shortID}g"
    print("@@@",gparam_specific_dir)
    print("####", subfoldername)

    # cr_dir(gparam_specific_dir)
    if os.path.exists(gparam_specific_dir):
        pass
    else:
        os.mkdir(gparam_specific_dir)

    clustexplicit = f"{gparam_specific_dir}/{os.path.basename(gparam_specific_dir)}.Explicit.cluster.info"
    clustshort = f"{gparam_specific_dir}/{os.path.basename(gparam_specific_dir)}.Short.cluster.info"
    clustoverlap = f"{gparam_specific_dir}/{os.path.basename(gparam_specific_dir)}.Overlap.cluster.info"
    clustdict = f"{gparam_specific_dir}/{os.path.basename(gparam_specific_dir)}.Dict.cluster.dict"
    
    print(clustexplicit)
    print(clustshort)
    print(clustoverlap)
    print(clustdict)

    with open(clustexplicit, "w") as f1, open(clustshort, "w") as f2, open(clustoverlap, "w") as f3, open(clustdict, "w") as f4:
        matrix_list = [mi for mi in os.listdir(subfolder) if mi.endswith(".matrix")]

        print(f"- Processing data from: {subfolder}")
        print(f"-- Total Number of input matrices: {len(matrix_list)}")
        process_matrices(matrix_list, f1, f2, f3, f4, subfolder, g_param_value, ExpandControl_param)

        print("# Output dictionary with cluster information: ", clustdict)

finish = time.perf_counter()
print(f"Finished in {round(finish-start, 2)} second(s)") # the result is 0.0 seconds just because our code have run concurrently, without waiting previous jobs to be completed
