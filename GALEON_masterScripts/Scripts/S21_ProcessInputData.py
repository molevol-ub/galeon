import os, shutil, subprocess, ast, sys
import pandas as pd
from collections import Counter
from operator import attrgetter

''' Parameteres and PATHs '''

GFF_dir = sys.argv[1] # "GFFs"
bedtoolsPATH = sys.argv[2] # PATH to bedtools
Feat2Grep = sys.argv[3] # Freature2Grep

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
				if ifile.endswith("_fam.gff3"):
					temp2.append(ifile)

				elif ifile.endswith("_fam.bed1"):
					temp2.append(ifile)

				elif ifile.endswith("_fam.bed2"):
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
	def __init__(self, scaffold, start, end, attribute, clusterID):
		self.scaffold = scaffold

		self.start = start
		self.end = end

		self.attribute = attribute
		self.clusterID = clusterID

	def __repr__(self):
		return repr([self.scaffold, self.start, self.end, self.attribute, self.clusterID])

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

			if colnum != 9:
				emsg = f"If this file: {i_file} has this format '{fileformat}', it should contain 9 columns separated by tabular delimiter. This file has {colnum} columns"
				print(emsg)
				raise ValueError("1")
				raise ValueError(emsg)
			else:
				return True

		elif fileformat == "bed2":
			DF = pd.read_csv(i_file, sep="\t")
			temp = DF.shape
			colnum = temp[1]

			if colnum != 4:
				emsg = f"If this file: {i_file} has this format '{fileformat}', it should contain 4 columns separated by tabular delimiter. This file has {colnum} columns"
				print(emsg)
				raise ValueError(emsg)
			else:
				return True

		elif fileformat == "bed1":
			DF = pd.read_csv(i_file, sep="\t")
			temp = DF.shape
			colnum = temp[1]

			if colnum != 3:
				emsg = f"If this file: {i_file} has this format '{fileformat}', it should contain 3 columns separated by tabular delimiter. This file has {colnum} columns"
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
				# print("Bedtools version: ", temp.stdout)
				pass

		else:
			emsg = f"The binary file '{i_binarypath}' does not exist or is not executable."
			raise ValueError(emsg)

	else:
		print("Bedtools already in path")
		print("Bedtools version: ", x.stdout)

CheckBedtools(bedtoolsPATH)

# F4: Grep 'gene' entries
def grep_entries(i_file, o_file, i_fileformat, i_featname=None):

	if i_fileformat == "gff3":
		grep_cmd = f"grep -P '\\t{i_featname}\\t' {i_file} > {o_file}"
		sortcluster_cmd = f"grep -P '\\t{i_featname}\\t' {i_file} | cut -f1,4,5,9 | sort -k1,1 -k2,2n | bedtools cluster -i - > {o_file}.bedclust\n"
		grep_count = f"grep -cP '\\t{i_featname}\\t' {o_file}"

		subprocess.run(grep_cmd, shell=True, text=True)
		subprocess.run(sortcluster_cmd, shell=True, text=True)
		gene_count = subprocess.run(grep_count, shell=True, text=True, capture_output=True).stdout.strip()

		o_msg = f"- Gene entries in {i_file}: {gene_count}"
		print(o_msg)

		return f"{o_file}.bedclust"

	elif i_fileformat == "bed2":
		catsort_cmd = f"cat {i_file} | sort -k1,1 -k2,2n | bedtools cluster -i - > {o_file}.bedclust\n"

		grep_count = f"wc -l {o_file}"

		subprocess.run(catsort_cmd, shell=True, text=True)
		gene_count = subprocess.run(grep_count, shell=True, text=True, capture_output=True).stdout.strip()

		o_msg = f"- Gene entries in {i_file}: {gene_count}"
		print(o_msg)

		return f"{o_file}.bedclust"

	elif i_fileformat == "bed1":
		catsort_cmd = f"cat {i_file} | sort -k1,1 -k2,2n | bedtools cluster -i - > {o_file}.bedclust\n"

		grep_count = f"cat {i_file} | wc -l"

		subprocess.run(catsort_cmd, shell=True, text=True)
		gene_count = subprocess.run(grep_count, shell=True, text=True, capture_output=True).stdout.strip()

		o_msg = f"- Gene entries in {i_file}: {gene_count}"
		print(o_msg)

		return f"{o_file}.bedclust"


#F5: Sort the grepped entries by Scaffold name and Start coordinate (col1, col4 (as integer); Check if start_coord > end_coord)
def load_and_sort(i_file, o_gff, i_genefamid):
	print("\n- Sorting the retrieving 'gene' entries by Scaffold name and Start coordinate")
	print(" - Input temp file: ", i_file)
	print(" - Output temp file: ", o_gff)

	loaded_gfflines = []

	with open(i_file, "r") as f1:

		for line in f1:
			line_s = line.strip().split("\t")
			scf_x, coord1, coord2, attr_x, cluster_ID = line_s[0], int(line_s[1]), int(line_s[2]), line_s[3], line_s[4]

			info = [scf_x, coord1, coord2, attr_x, cluster_ID]

			if coord1 >= coord2:
				print("!!! Error (start_coord >= end_coord), it will be fixed by reversing the order (for ex: 10 5 => 5 10). The user must check if that solution makes sense for the info stored in this GFF.")
				print("!!!", line_s)
				loaded_gfflines.append(LoadBED(scf_x, coord2, coord1, attr_x, cluster_ID))
			else:
				loaded_gfflines.append(LoadBED(scf_x, coord1, coord2, attr_x, cluster_ID))

	print(f"\n- Looking for overlapping genes in {o_gff}")
	
	# Create en empty dictionary that will hold all the cluster IDs
	gene_dict = {g.scaffold : {"Clusters" : []} for g in loaded_gfflines}

	# Fill it with cluster IDs list
	for g in loaded_gfflines:
		gene_dict[g.scaffold]["Clusters"].append(g.clusterID)
	
	# Check whether some cluster IDs are present more than once => there are several genes (at least 2) that are overlapping
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
		print(f"-- This file has {len(temp_dict_2)} overlapping sets of genes --")
		outgene_entries = []
		D_bedclust2genename = {}

		for iii in loaded_gfflines:
			if iii.clusterID not in temp_dict_2.keys():
				outgene_entries.append(iii)

			else:
				temp_dict_2[iii.clusterID]["info"].append(iii)

		### checkpoint
		for k, v in temp_dict_2.items():
			if v["gene_count"] != len(v["info"]):
				msg = "Something went wrong, probably with BEDTOOLS cluster"
				raise ValueError(msg)
		##

		control_counter = len(loaded_gfflines)
		###	

		for k, v in temp_dict_2.items():
			print("----> Gene_Overlap case #{} : {} genes".format(k, v["gene_count"]))
			for cigene in v["info"]:
				print("\t", cigene)

			t_value = v["gene_count"] -1
			control_counter -= t_value

			# this will be used as a out cluster name # GG{newID}_{overlapcase}_{gene_count}
			newID = max(overlapcases_numset) + 1
			overlapcases_numset.add(newID)

			bedgroup_name = "GG#{}_{}_{}".format(newID, k, v["gene_count"])
			D_bedclust2genename[bedgroup_name] = v["info"]
			print("----> This set of collapsed genes will be named as: ", bedgroup_name, "\n")

			templist_coord1 = []
			templist_coord2 = []

			for gene_entry in v["info"]:
				templist_coord1.append(gene_entry.start)
				templist_coord2.append(gene_entry.end)

			bedgroup_start, bedgroup_end = min(templist_coord1), max(templist_coord2)

			outgene_entries.append(LoadBED(gene_entry.scaffold, bedgroup_start, bedgroup_end, bedgroup_name, k))

			with open(f"{GFF_dir}/{i_genefamid}_OverlappingGenes.txt", "w") as o1, open(f"{GFF_dir}/{i_genefamid}_OverlappingGenes.dict", "w") as o2:
				header = "\t".join(["#Scaffold", "start", "end", "attribute", "clusterID"])
				print(header, end="\n", file=o1)

				for k,v in temp_dict_2.items():
					for vi in v["info"]:
						if "ID=" in vi.attribute:
							vi.attribute = vi.attribute.split(";")[0].replace("ID=","")
							print(vi.scaffold, vi.start, vi.end, vi.attribute, vi.clusterID, sep="\t", file=o1)
						else:
							print(vi.scaffold, vi.start, vi.end, vi.attribute, vi.clusterID, sep="\t", file=o1)
				o2.write(str(D_bedclust2genename))
			
		if len(outgene_entries) == control_counter:
			print("# INPUT lines from gff temp file: ", len(loaded_gfflines))
			print("# OUTPUT lines: ", len(outgene_entries))
		else:
			msg = "# Something is wrong"
			raise ValueError(msg)

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


# STEP 2: Check that every gene entry in the GFF file has the following format (start_coordinate > end_coordinate)
# Afterwards, cluster overlapping coordinates and generate a *temp file that will next be processed

msg = "# | S2.1.3 Step | Checking Annotation Files"
print(msg)

out_tempfiles = []

if annotfiles_format == "gff3":
	for i in annotfiles_list:
		GeneFamID = i.split("_fam.")[0] + "_fam"
		coord_file = f"{GFF_dir}/{i}" # GFF
		CheckAnnotFile(coord_file)
			
		temp_file = f"{i}.temp" # temp file

		# remove these temporary files (if they are already present)
		check_dir(GFF_dir, "file", temp_file)

		temp_file = f"{GFF_dir}/{temp_file}"

		# Grep genes, sort and cluster
		temp_file = grep_entries(coord_file, temp_file, annotfiles_format, Feat2Grep)
		out_tempfiles.append(temp_file)


	for tempfile in out_tempfiles:
		temp = os.path.basename(tempfile)
		GeneFamID = temp.split("_fam.")[0] + "_fam"
		print(tempfile)
		collapsed_tempfile = ".".join(tempfile.split(".")[:-1]) + ".collapsed.temp"
		print(collapsed_tempfile)
		load_and_sort(tempfile, collapsed_tempfile, GeneFamID)

elif annotfiles_format == "bed2":
	for i in annotfiles_list:
		GeneFamID = i.split("_fam.")[0] + "_fam"

		coord_file = f"{GFF_dir}/{i}" # GFF
		CheckAnnotFile(coord_file)
			
		temp_file = f"{i}.temp" # temp file

		# remove these temporary files (if they are already present)
		check_dir(GFF_dir, "file", temp_file)

		temp_file = f"{GFF_dir}/{temp_file}"

		temp_file = grep_entries(coord_file, temp_file, annotfiles_format)
		out_tempfiles.append(temp_file)		

	for tempfile in out_tempfiles:
		temp = os.path.basename(tempfile)
		GeneFamID = temp.split("_fam.")[0] + "_fam"
		print(tempfile)
		collapsed_tempfile = ".".join(tempfile.split(".")[:-1]) + ".collapsed.temp"
		print(collapsed_tempfile)
		load_and_sort(tempfile, collapsed_tempfile, GeneFamID)

elif annotfiles_format == "bed1":
	for i in annotfiles_list:
		GeneFamID = i.split("_fam.")[0] + "_fam"
		GeneFam_name = i.replace("_fam.bed1","")

		coord_file = f"{GFF_dir}/{i}" # GFF
		CheckAnnotFile(coord_file)
			
		temp_file = f"{i}.temp" # temp file

		# remove these temporary files (if they are already present)
		check_dir(GFF_dir, "file", temp_file)

		# Add arbitrary genenames
		input_bedfile = f"{GFF_dir}/{i}"
		itemp = i + ".genenames"
		output_bedfile = f"{GFF_dir}/{itemp}"
		temp2_file = f"{GFF_dir}/{itemp}.temp"

		with open(input_bedfile) as iifile, open(output_bedfile, "w") as ofile:
			for inum, iline in enumerate(iifile):
				line_s = iline.strip().split("\t")
				Arbitrary_GeneID = f"{GeneFam_name}_{inum}"
				line_s += [Arbitrary_GeneID]
				jline = "\t".join(line_s) + "\n"
				ofile.write(jline)

		temp2_file = grep_entries(output_bedfile, temp2_file, annotfiles_format)
		out_tempfiles.append(temp2_file)		

	for tempfile in out_tempfiles:
		temp = os.path.basename(tempfile)
		GeneFamID = temp.split("_fam.")[0] + "_fam"
		print(tempfile)
		collapsed_tempfile = ".".join(tempfile.split(".")[:-1]) + ".collapsed.temp"
		print(collapsed_tempfile)
		load_and_sort(tempfile, collapsed_tempfile, GeneFamID)

else:
	err_msg = f"Unknown input annotation file format: {annotfiles_format}"
	raise ValueError(err_msg)
