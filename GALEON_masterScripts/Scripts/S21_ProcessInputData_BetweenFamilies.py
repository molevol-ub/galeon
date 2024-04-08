import os, shutil, subprocess, ast, sys
import pandas as pd
from collections import Counter
from operator import attrgetter

''' Parameteres and PATHs '''

GFF_dir = sys.argv[1] # "GFFs"
GFF_mergedir = f"{GFF_dir}/merged_dir"
bedtoolsPATH = sys.argv[2] # PATH to bedtools

# bash_cmds = "sort_merge.sh"

overlapcases_numset = set([0])

''' Functions '''
# F1: create or remove directories in working directory

def cr_dir(i_dir): 
    # Fx to create/overwrite a directory
    
    if i_dir in os.listdir():
        # print(f"Removing old directory and creating a new one: '{i_dir}'")
        shutil.rmtree(i_dir)
        os.mkdir(i_dir)
    
    else:
        # print(f"Creating a new directory: '{i_dir}'")
        os.mkdir(i_dir)		

# F: check the format of the input data
def check_input_format(i_list):
	format_set = set()

	for ii in i_list:
		temp = ii.split(".")[-1]
		format_set.add(temp)

	if len(format_set) == 0:
		raise ValueError("The input Annotation files must end in *gff3 / *bed1 / *bed2 ")
	
	elif len(format_set) != 1:
		raise ValueError("All the input Annotation files must have the same format: *gff3 / *bed1 / *bed2 ")
	else:
		return list(format_set)

# F2: check if a given folder is present in working directory
def check_dir(i_dirname, i_opt, i_file=""): 
	wkd_dirlist = [i[0] for i in os.walk(".")]

	if i_opt == "check":
		print(f"# Checking if the folder: '{i_dirname}' is present in the current directory: '.'")
		
		if i_dirname in wkd_dirlist:
			# print("## - Folder found! It will be overwritten")
			cr_dir(i_dirname) #, "remove")
			# print("## --- Folder removed")

			# cr_dir(temp) #, "create")

		else:
			# print("## - Folder not found! It will be created")
			cr_dir(i_dirname) #, "create")
			# print("## --- Folder created")

	elif i_opt == "content":
		
		content_list = []

		try:
			# temp2 = [ifile for ifile in os.listdir(i_dirname) if ifile.endswith(".gff3")]
			temp2 = []

			for ifile in os.listdir(i_dirname):
				if ifile.endswith(".gff3"):
					temp2.append(ifile)

				elif ifile.endswith(".bed1"):
					temp2.append(ifile)

				elif ifile.endswith(".bed2"):
					temp2.append(ifile)

				else:
					pass


			if len(temp2) == 0:
				if len(os.listdir(i_dirname)) == 0:
					raise ValueError(f"! Error: the input directory is empty", i_dirname)
				else:
					raise ValueError(f"! Error. Check that your GFF files are named as (*_fam.gff3 or *_fam.bed1 or *_fam.bed2)")

			else:
				input_data_format = check_input_format(temp2)[0]

				print(f"The input files from the a annotation directory '{i_dirname}' were provided in this format '{input_data_format}'")

				for afile in temp2:
		
					if afile.endswith("_fam.gff3"): 
						content_list.append(afile)
					
					elif afile.endswith("_fam.bed1"): 
						content_list.append(afile)

					elif afile.endswith("_fam.bed2"): 
						content_list.append(afile)
					else:
						pass

				if len(content_list) == 0:
					raise ValueError(f"! Error. Check that your GFF files are named as (*_fam.gff3 or *_fam.bed1 or *_fam.bed2)")
		
		except FileNotFoundError as e:
			msg = f"! Error. Input directory '{i_dirname}' not found"
			raise FileNotFoundError(msg)

		return content_list, input_data_format

	elif i_opt == "file":
		rmfile = f"{i_dirname}/{i_file}"
		if i_file in os.listdir(i_dirname):
			os.remove(rmfile)
		else:
			pass


# F3: Class to load the gff lines and assign a description to each element
class LoadGFF:
	def __init__(self, scaffold, source, feature, start, end, score, strand, frame, attribute):
		self.scaffold = scaffold
		self.source = source 
		self.feature = feature

		self.start = start
		self.end = end
		self.score = score

		self.strand = strand
		self.frame = frame
		self.attribute = attribute

	def __repr__(self):
		return repr([self.scaffold, self.source, self.feature, self.start, self.end, self.score, self.strand, self.frame, self.attribute])

class LoadBED:
	def __init__(self, scaffold, start, end, attribute, famID, clusterID):
		self.scaffold = scaffold

		self.start = start
		self.end = end

		self.attribute = attribute
		self.famID = famID
		self.clusterID = clusterID

	def __repr__(self):
		return repr([self.scaffold, self.start, self.end, self.attribute, self.famID, self.clusterID])

# Fx to check input file format
def CheckAnnotFile(i_file):
	print(i_file)
	fileformat = i_file.split(".")[-1]

	if fileformat not in ["gff3", "bed1", "bed2"]:
		emsg = f"Unknown format '{fileformat}' for this file: {i_file}"
		raise ValueError(emsg)

	else:
		if fileformat == "gff3":
			DF = pd.read_csv(i_file, sep="\t")
			temp = DF.shape
			colnum = temp[1]

			if colnum != 9+1:
				emsg = f"If this file: {i_file} has this format '{fileformat}', it should contain 9 columns separated by tabular delimiter (+1 FamilyID column when comparing two gene families). This file has {colnum} columns"
				print(emsg)
				raise ValueError("1")
				raise ValueError(emsg)
			else:
				return True

		elif fileformat == "bed2":
			DF = pd.read_csv(i_file, sep="\t")
			temp = DF.shape
			colnum = temp[1]

			if colnum != 4+1:
				emsg = f"If this file: {i_file} has this format '{fileformat}', it should contain 4 columns separated by tabular delimiter (+1 FamilyID column when comparing two gene families). This file has {colnum} columns"
				print(emsg)
				raise ValueError(emsg)
			else:
				return True

		elif fileformat == "bed1":
			DF = pd.read_csv(i_file, sep="\t")
			temp = DF.shape
			colnum = temp[1]

			if colnum != 4+1:
				emsg = f"If this file: {i_file} has this format '{fileformat}', it should contain 3 columns separated by tabular delimiter (+1 FamilyID column when comparing two gene families). This file has {colnum} columns"
				print(emsg)
				raise ValueError(emsg)
			else:
				return True

# Fx to check that bedtools is available
def CheckBedtools(i_binarypath=""):
	# If the bedtools is already in path, nothing is done
	x = subprocess.run("bedtools --version", shell=True, text=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE)

	if x.stdout == "":
		# Check if the binary file exists and if it is executable
		if os.path.exists(i_binarypath) and os.access(i_binarypath, os.X_OK):
			# Get the current PATH and append the binary's directory to it
			os.environ['PATH'] = f"{os.path.dirname(i_binarypath)}:{os.environ['PATH']}"

			# Now you should be able to run the binary using its name directly
			temp = subprocess.run("bedtools --version", shell=True, text=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
			print("User specified PATH to bedtools: ", i_binarypath)
			print("Bedtools version: ", temp.stdout)

			if temp.stdout == "":
				emsg = f"The bedtools binary file does not exist in this directory '{i_binarypath}' or is not executable"
				raise ValueError(emsg)
			else:
				print("Bedtools version: ", temp.stdout)

		else:
			emsg = f"The binary file '{i_binarypath}' does not exist or is not executable."
			raise ValueError(emsg)

	else:
		print("Bedtools already in path")
		print("Bedtools version: ", x.stdout)
CheckBedtools(bedtoolsPATH)

# F4: Grep 'gene' entries
def Grep_and_BedCluster(i_file, o_file, i_format):
	
	if i_format == "gff3":
		# grep_cmd = f"grep -P '\\tgene\\t' {i_file} > {o_file}"
		sortcluster_cmd = f"grep -P '\\tgene\\t' {i_file} | cut -f1,4,5,9,10 | sort -k1,1 -k2,2n | {bedtoolsPATH} cluster -i - > {o_file}.bedclust\n"
		grep_count = f"grep -cP '\\tgene\\t' {o_file}"

		# subprocess.run(grep_cmd, shell=True, text=True)
		subprocess.run(sortcluster_cmd, shell=True, text=True)
		gene_count = subprocess.run(grep_count, shell=True, text=True, capture_output=True).stdout.strip()

		o_msg = f"- Gene entries in {i_file}: {gene_count}"
		print(o_msg)

		return f"{o_file}.bedclust"

	elif i_format == "bed1":
		outcmd = f"cat {i_file} | sort -k1,1 -k2,2n | {bedtoolsPATH} cluster -i - > {o_file}.bedclust\n"
		# print(outcmd, o_file)
		grep_count = f"cat {i_file} | wc -l"

		subprocess.run(outcmd, shell=True, text=True)
		gene_count = subprocess.run(grep_count, shell=True, text=True, capture_output=True).stdout.strip()

		# o_msg = f"- Gene entries in {i_file}: {gene_count}"
		o_msg = f"- Gene entries in merged file: {gene_count}"
		# print(o_msg)

		return f"{o_file}.bedclust"


	elif i_format == "bed2":
		catsort_cmd = f"cat {i_file} | sort -k1,1 -k2,2n | {bedtoolsPATH} cluster -i - > {o_file}.bedclust\n"
		# print(catsort_cmd)
		grep_count = f"cat {i_file} | wc -l"

		subprocess.run(catsort_cmd, shell=True, text=True)
		gene_count = subprocess.run(grep_count, shell=True, text=True, capture_output=True).stdout.strip()

		# o_msg = f"- Gene entries in {i_file}: {gene_count}"
		o_msg = f"- Gene entries in merged file: {gene_count}"
		print(o_msg)

		return f"{o_file}.bedclust"


#F5: Sort the grepped entries by Scaffold name and Start coordinate (col1, col4 (as integer); Check if start_coord > end_coord)
def load_and_sort(i_file, o_gff):
	print("\n- Sorting the retrieved 'gene' entries by 'Scaffold name' and 'Start coordinate'")
	print(" - Input temp file: ", i_file)
	print(" - Output temp file: ", o_gff)

	loaded_gfflines = []

	with open(i_file, "r") as f1:

		for line in f1:
			line_s = line.strip().split("\t")
			scf_x, coord1, coord2, attr_x, GeneFam_ID, cluster_ID = [*line_s]
			coord1, coord2 = int(coord1), int(coord2)

			# info = [scf_x, coord1, coord2, attr_x, GeneFam_ID, cluster_ID]

			if coord1 >= coord2:
				print("!!! Error (start_coord >= end_coord), it will be fixed by reversing the order (for ex: 10 5 => 5 10). The user must check if that solution makes sense for the info stored in this GFF.")
				print("!!!", line_s)
				loaded_gfflines.append(LoadBED(scf_x, coord2, coord1, attr_x, GeneFam_ID, cluster_ID))
			else:
				loaded_gfflines.append(LoadBED(scf_x, coord1, coord2, attr_x, GeneFam_ID, cluster_ID))

	# print(f"\n- Looking for overlapping genes in {o_gff}")
	print(f"\n- Looking for overlapping genes")
	
	# Create en empty dictionary that will hold all the cluster IDs
	gene_dict = {g.scaffold : {"Clusters" : []} for g in loaded_gfflines}

	# Fill it with cluster IDs list
	for g in loaded_gfflines:
		gene_dict[g.scaffold]["Clusters"].append(g.clusterID)
	
	# Check whether some cluster IDs are present more than once => there are several genes (at least 2) that are overlapped
	temp_dict_2 = {}

	for k_scf, v_values in gene_dict.items():
		set_length = len(set(v_values["Clusters"]))
		list_length = len(v_values["Clusters"])
		
		if set_length == list_length:
			pass
		else:

			temp_dict = dict(Counter(v_values["Clusters"]))

			for k,v in temp_dict.items():
				if v != 1:
					temp_dict_2[k] = {"gene_count" : v, "info" : []} 

			# print(temp_dict_2)
	
	if len(temp_dict_2) == 0:

		print(" - No overlapping genes have been found")

		with open(o_gff, "w") as o1:

			for iline in loaded_gfflines:
				x = ast.literal_eval(str(iline))
				newx = []

				for xx in x: newx.append(str(xx))

				jline = "\t".join(newx) + "\n"
				o1.write(jline)
			# print("## Done here ##")

	else:
		D_bedclust2genename = {}

		print(f"-- This file has {len(temp_dict_2)} overlapping sets of genes --")
		print(f"--- Check this 'Cluster IDs' in *temp files: {list(temp_dict_2.keys())}")
		outgene_entries = []

		for iii in loaded_gfflines:
			if iii.clusterID not in temp_dict_2.keys():
				outgene_entries.append(iii)

			else:
				temp_dict_2[iii.clusterID]["info"].append(iii)

		### checkpoint
		for k, v in temp_dict_2.items():
			if v["gene_count"] != len(v["info"]):
				msg = "something went wrong, probably with BEDTOOLS cluster"
				raise ValueError(msg)
		##
		print(temp_dict_2)

		control_counter = len(loaded_gfflines)
		###	
		for k, v in temp_dict_2.items():
			# print(k, "@@@@")
			# print(v) # @@@@
			print("---- Gene_Overlap case #{} : {} genes".format(k, v["gene_count"]))
			for cigene in v["info"]:
				print("\t", cigene)

			t_value = v["gene_count"] -1
			control_counter -= t_value

			# this will be used as a out cluster name # GG{newID}_{overlapcase}_{gene_count}
			newID = max(overlapcases_numset) + 1
			overlapcases_numset.add(newID)
            
            # Create a label (hint) to know which families were merged
			mlabel_temp = []
			for iinfo in v["info"]:
				if iinfo.famID not in mlabel_temp:
					mlabel_temp.append(iinfo.famID)
				else:
					pass
			mlabel = "--".join(mlabel_temp)

			bedgroup_name = "{}_GG#{}_{}_{}".format(mlabel, newID, k, v["gene_count"])
			D_bedclust2genename[bedgroup_name] = v["info"]
			print("----> This set of collapsed genes will be named as: ", bedgroup_name, "\n")

			templist_coord1 = []
			templist_coord2 = []

			for gene_entry in v["info"]:
				templist_coord1.append(gene_entry.start)
				templist_coord2.append(gene_entry.end)

			bedgroup_start, bedgroup_end = min(templist_coord1), max(templist_coord2)

			outgene_entries.append(LoadBED(gene_entry.scaffold, bedgroup_start, bedgroup_end, bedgroup_name, gene_entry.famID, k))


		if len(outgene_entries) == control_counter:
			print("# INPUT lines from gff temp file: ", len(loaded_gfflines))
			print("# OUTPUT lines: ", len(outgene_entries))
		else:
			msg = "# Something is wrong"
			raise ValueError(msg)

		with open(f"{GFF_mergedir}/merged_fam_OverlappingGenes.dict", "w") as o2:
			o2.write(str(D_bedclust2genename))

		# Exporting results
		outgene_entries.sort(key=attrgetter("scaffold", "start"))

		with open(o_gff, "w") as o1:
			for iline in outgene_entries:
				x = ast.literal_eval(str(iline))
				newx = []

				for xx in x: newx.append(str(xx))

				jline = "\t".join(newx) + "\n"
				o1.write(jline)

			print("## Done here ##")


''' Processing '''
# STEP1: Check the content of 'GFFs' folder

msg = "# | S2.1.1 Step | Checking Annotation Directory"
print(msg)

annotfiles_list, annotfiles_format = check_dir(GFF_dir, "content") # this folder must contain coord_files to be analyzed

msg = "# | S2.1.2 Step | The following annotation files will be analyzed"
print(msg)

for i in annotfiles_list:
	print(f"- {i}")

print("\n")


# Merge GFF files
if os.path.exists(GFF_mergedir) == True:
	shutil.rmtree(GFF_mergedir)
	os.mkdir(GFF_mergedir)
else:
	os.mkdir(GFF_mergedir)

# Load each gff as pandas dataframe

def load_gff_as_pd(i_file, i_dir, i_format):
	if i_format == "gff3":
		Df_out = pd.read_csv(f"{i_dir}/{i_file}", sep="\t", header=None)
		
		FamName = i_file.split("_")[0] # Example: 'GR_fam.gff3' => 'GR'
		Df_out[9] = FamName

		return FamName, Df_out

	elif i_format == "bed2":
		Df_out = pd.read_csv(f"{i_dir}/{i_file}", sep="\t", header=None)

		FamName = i_file.split("_")[0] # Example: 'GR_fam.bed2' => 'GR'
		Df_out[4] = FamName
		
		return FamName, Df_out

	elif i_format == "bed1":
		Df_out = pd.read_csv(f"{i_dir}/{i_file}", sep="\t", header=None)
		
		FamName = i_file.split("_")[0] # Example: 'GR_fam.bed1' => 'GR'
		ArbitGeneNames = []

		for idxx in range(0,len(Df_out)):
			AGeneID = f"{FamName}_{idxx}"
			ArbitGeneNames.append(AGeneID)

		Df_out[4] = ArbitGeneNames
		Df_out[5] = FamName

		return FamName, Df_out

# Merge two gff
# - One extra field with gene family name is going to be added
def merge_CoordFiles(i_fileslist, i_dir, o_dir, i_format):
	if i_format == "gff3":
		merged_df = None
		
		GeneFam_list = {}
		
		for file in i_fileslist:
			famname, temp = load_gff_as_pd(file, i_dir, i_format)
			GeneFam_list[famname] = temp
		
		temp = len(GeneFam_list.values())
		if temp != 2:
			msg = f"! ERROR: Only two families are allowed to be as input for merging! The use provided: {temp}"
			print(msg)
			raise ValueError(msg)

		else:
			merged_df = pd.concat(GeneFam_list.values())
			merged_df.to_csv(f"{o_dir}/merged_fam.gff3", sep = "\t", index=False, header=None)

	
	elif i_format == "bed2":
		merged_df = None
		
		GeneFam_list = {}
		
		for file in i_fileslist:
			famname, temp = load_gff_as_pd(file, i_dir, i_format)
			GeneFam_list[famname] = temp


		temp = len(GeneFam_list.values())
		if temp != 2:
			msg = f"! ERROR: Only two families are allowed to be as input for merging! The use provided: {temp}"
			print(msg)
			raise ValueError(msg)

		else:
			merged_df = pd.concat(GeneFam_list.values())
			merged_df.to_csv(f"{o_dir}/merged_fam.bed2", sep = "\t", index=False, header=None)

	elif i_format == "bed1":
		merged_df = None
		
		GeneFam_list = {}
		
		for file in i_fileslist:
			famname, temp = load_gff_as_pd(file, i_dir, i_format)
			GeneFam_list[famname] = temp


		temp = len(GeneFam_list.values())
		if temp != 2:
			msg = f"! ERROR: Only two families are allowed to be as input for merging! The use provided: {temp}"
			print(msg)
			raise ValueError(msg)

		else:
			merged_df = pd.concat(GeneFam_list.values())
			merged_df.to_csv(f"{o_dir}/merged_fam.bed1", sep = "\t", index=False, header=None)

merge_CoordFiles(annotfiles_list, GFF_dir, GFF_mergedir, annotfiles_format)

annotfiles_list = [i for i in os.listdir(GFF_mergedir) if i.endswith(f".{annotfiles_format}")]

# STEP 2: Check that every gene entry in the GFF file has the following format (start_coordinate > end_coordinate)
# Afterwards, cluster overlapping coordinates and generate a *temp file that will next be processed

msg = "# | S2.1.3 Step | Checking Annotation Files"
print(msg)

out_tempfiles = []

if annotfiles_format == "gff3":
	for i in annotfiles_list:
		coord_file = f"{GFF_mergedir}/{i}" # GFF
		CheckAnnotFile(coord_file)
		temp_file = f"{i}.temp" # temp file

		# remove these temporary files (if they are already present)
		check_dir(GFF_dir, "file", temp_file)

		temp_file = f"{GFF_mergedir}/{temp_file}"
		# out_gff = f"{GFF_sorted_dir}/{out_gff_temp}"

		temp_file = Grep_and_BedCluster(coord_file, temp_file, annotfiles_format)
	# 	# bashcmd_file.write(bashcmd)
		# print(temp_file)
		out_tempfiles.append(temp_file)

	# Sort and cluster genes
	for tempfile in out_tempfiles:

		collapsed_tempfile = ".".join(tempfile.split(".")[:-1]) + ".collapsed.temp"

		load_and_sort(tempfile, collapsed_tempfile)

elif annotfiles_format == "bed2":
	for i in annotfiles_list:
		print(i)
		coord_file = f"{GFF_mergedir}/{i}" # GFF
		CheckAnnotFile(coord_file)
		temp_file = f"{i}.temp" # temp file

		# remove these temporary files (if they are already present)
		check_dir(GFF_dir, "file", temp_file)

		temp_file = f"{GFF_mergedir}/{temp_file}"
		# out_gff = f"{GFF_sorted_dir}/{out_gff_temp}"

		temp_file = Grep_and_BedCluster(coord_file, temp_file, annotfiles_format)
	# 	# bashcmd_file.write(bashcmd)
		# print(temp_file)
		out_tempfiles.append(temp_file)

	# Sort and cluster genes
	for tempfile in out_tempfiles:

		collapsed_tempfile = ".".join(tempfile.split(".")[:-1]) + ".collapsed.temp"

		load_and_sort(tempfile, collapsed_tempfile)

elif annotfiles_format == "bed1":
	for i in annotfiles_list:
		print(i)
		coord_file = f"{GFF_mergedir}/{i}" # GFF
		CheckAnnotFile(coord_file)
		temp_file = f"{i}.temp" # temp file

		# remove these temporary files (if they are already present)
		check_dir(GFF_dir, "file", temp_file)

		temp_file = f"{GFF_mergedir}/{temp_file}"
		# out_gff = f"{GFF_sorted_dir}/{out_gff_temp}"

		temp_file = Grep_and_BedCluster(coord_file, temp_file, annotfiles_format)
	# 	# bashcmd_file.write(bashcmd)
		# print(temp_file)
		out_tempfiles.append(temp_file)

	# Sort and cluster genes
	for tempfile in out_tempfiles:
		print("!!!!!!!!",tempfile)
		collapsed_tempfile = ".".join(tempfile.split(".")[:-1]) + ".collapsed.temp"

		load_and_sort(tempfile, collapsed_tempfile)
