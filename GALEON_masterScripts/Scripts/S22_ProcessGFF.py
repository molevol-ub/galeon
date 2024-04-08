import os, sys, subprocess, shutil
import pandas as pd
from itertools import combinations_with_replacement
import time
start = time.perf_counter()


''' Input '''

inputdir = sys.argv[1]
outputdir = sys.argv[2]

temp = os.path.abspath(inputdir) # wkdirectory
wkd = os.path.dirname(temp) # wkdirectory


# ''' Parameteres '''

# ''' Functions '''
        
def cr_dir(i_dir): 
    # Fx to create/overwrite a directory
    
    if os.path.isdir(i_dir):
        # print(f"Removing old directory and creating a new one: '{i_dir}'")
        shutil.rmtree(i_dir)
        os.mkdir(i_dir)
    
    else:
        # print(f"Creating a new directory: '{i_dir}'")
        os.mkdir(i_dir)

def Process_input(i_file, i_dir, i_opt):
    x_dict = {}
    x_list = []
    
    # cut scaffold field
    cut_cmd = f"cut -f1 {i_dir}/{i_file} | sort | uniq"
    af_scaffolds = subprocess.run(cut_cmd, shell=True, text=True, capture_output=True).stdout.strip().split("\n")
    af_gene_dict = {p1i : {} for p1i in af_scaffolds}
    
    if i_opt == "GFF":
        scaf_loc, field_type, start_loc, end_loc, orientation, atrib_loc, newcol_loc = 0, 2, 3, 4, 6, 8, 9

        # import gff information
        with open(f"{i_dir}/{i_file}", "r") as ifile:
            for line in ifile:
                line = line.strip().split("\t")

                if line[field_type] == "gene":
                    scfID = line[scaf_loc]
                    gstart = int(line[start_loc])
                    gend = int(line[end_loc])
                    gname = atrib_to_genename(line[atrib_loc]) # extract gene name from GFF information
                    
                    if gstart < gend:
                        af_gene_dict[scfID][gname] = [gstart, gend]
                    
                    elif gstart > gend:
                        
                        if line[orientation] == "-":
                            msg = f"- Flipping start/end coordinates ({gstart}/{gend}) for gene '{gname}' in scf '{scfID}'"
                            print(msg)
                            af_gene_dict[scfID][gname] = [gend, gstart]
                        
                        elif line[orientation] == "+":
                            warning = f"START coordinate {gstart} > END coordinate {gend} for gene '{gname}' in scf '{scfID}' with orientation '{line[orientation]}'!!! CHECK THAT !!!"
                            raise ValueError(warning)
                        else:
                            print(f"WARNING: Unknown orientation {line[orientation]}")
                    
                    elif gstart == gend:
                        warning = f"SAME START/END COORDINATES ({gstart},{gend}) for gene '{gname}' in scf '{scfID}'"
                        raise ValueError(warning)
                else:
                    pass
                
        return af_scaffolds, af_gene_dict

    if i_opt == "bed":
        scaf_loc, start_loc, end_loc, atrib_loc, bedcluster_loc, newcol_loc = 0, 1, 2, 3, 4, 5

        # import bed information

        with open(f"{i_dir}/{i_file}", "r") as ifile:
            for line in ifile:
                line = line.strip().split("\t")

                scfID = line[scaf_loc]
                gstart = int(line[start_loc])
                gend = int(line[end_loc])
                gname = atrib_to_genename(line[atrib_loc]) # extract gene name from GFF information
                    
                if gstart < gend:
                    af_gene_dict[scfID][gname] = [gstart, gend]
                
                elif gstart > gend:
                    msg = f"- Flipping start/end coordinates ({gstart}/{gend}) for gene '{gname}' in scf '{scfID}'"
                    print(msg)
                    af_gene_dict[scfID][gname] = [gend, gstart]
                        
                elif gstart == gend:
                    warning = f"SAME START/END COORDINATES ({gstart},{gend}) for gene '{gname}' in scf '{scfID}'"
                    raise ValueError(warning)
            else:
                pass
                
        return af_scaffolds, af_gene_dict

def atrib_to_genename(i_atrib):
    genename = i_atrib.split(";")[0].replace("ID=","")
    return genename

def combo_genes(i_genelist, i_opt):
    o_list = []
    
    if i_opt == "nr":
        o_list = list(combinations(i_genelist, 2))
    
    elif i_opt == "allcomb":
        o_list = list(combinations_with_replacement(i_genelist, 2))
        
    return o_list
    

def Compute_pairwise_PhysDist(i_dict, i_scf, i_afname, i_iodict, o_dir):

    # print(f"# --- Scf: {i_scf}")
    
    outname = f"{o_dir}/{i_iodict[i_afname]}/{i_iodict[i_afname]}.{i_scf}.matrix"
    
    i_genelist = list(i_dict[i_scf].keys()) # gene list from a given scaffold
    
    outheader = "\t".join(i_genelist) # matrix header 
    
    i_geneallcomb = combo_genes(i_genelist, "allcomb") # combinations of gene pairs
    
    # Compute pw distance and fill a dictionary

    pwdist_dict = {k1 : {} for k1 in i_genelist}

    for k,v in pwdist_dict.items():
        for kk in i_genelist:
            v[kk] = None

    for igcomb in i_geneallcomb:

        gene1 = igcomb[0]
        gene2 = igcomb[1]
            
        if gene1 == gene2:
            pwdist = [gene1,gene2,0]
            pwdist_dict[gene1][gene2] = 0
        
        else:
            
            #coordinates
            gene1_end = i_dict[i_scf][gene1][1] - 1  # convert coord from 1-based (gff) to 0-based (python) system
            gene2_start = i_dict[i_scf][gene2][0] - 1 # same here

            #print(i_dict[i_scf][gene1][1], i_dict[i_scf][gene2][0])
            mindist =  gene2_start - gene1_end + 1 # gene1 [0-5]; gene2 [10-20]; dist(gene2-gene1) = (5, 10] <=> [6,10] <=> 10 - 5 + 1

            if mindist < 1:
                print("# Case", igcomb)
                raise ValueError("# Distance < 1 !!!!")
            #info
            #info = f"{gene1}\t{gene2}\t{gene1_end}\t{gene2_start}\t{mindist}"

            pwdist_dict[gene1][gene2] = mindist
            pwdist_dict[gene2][gene1] = mindist

            if mindist < 0: 
                print("####")
                print(gene1, gene2, mindist)
                print(gene1, i_dict[i_scf][gene1])
                print(gene2, i_dict[i_scf][gene2])

    # list to dataframe
    out_df = pd.DataFrame.from_dict(pwdist_dict)

    del pwdist_dict
    
    if len(out_df) == 1:
        print(f"- Warning: Skipping scf {i_scf} from annotation file {i_afname}, only 1 gene has been found")
        pass
    else:
        out_df.to_csv(outname, sep="\t")


''' Processing '''

# Create output directory 
cr_dir(outputdir)
my_cwd = os.getcwd()
path_to_outdir = f"{my_cwd}/{outputdir}"

# Annotation files list (GFF, GTF, BED)
annotfiles_list = [i for i in os.listdir(inputdir) if i.endswith(".collapsed.temp")]
print(annotfiles_list)

IO_namesdict = {i : ".".join(i.split(".")[:-1]).replace(".collapsed","_matrices") for i in annotfiles_list}

print(IO_namesdict)

# Create subdirectories
for temp_filename, out_matrixdir in IO_namesdict.items():
    matrix_subfolder = f"{path_to_outdir}/{out_matrixdir}"
    os.mkdir(matrix_subfolder)



# Run the script over each file (and scaffold)

for afile in annotfiles_list:

    print(f"# Processing {afile} file")
    af_obj = Process_input(afile, inputdir, "bed")

    scflist, af_dict = af_obj[0], af_obj[1]
    for i_scf in scflist:
        Compute_pairwise_PhysDist(af_dict, i_scf, afile, IO_namesdict, outputdir)


finish = time.perf_counter()
print(f"Finished in {round(finish-start, 2)} second(s)") # the result is 0.0 seconds just because our code have run concurrently, without waiting previous jobs to be completed




