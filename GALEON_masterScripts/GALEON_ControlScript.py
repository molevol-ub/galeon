import argparse, ast, subprocess, os, shutil

temp = subprocess.run("which GALEON_ControlScript.py", shell=True, capture_output=True, text=True).stdout.strip()
PATHtoGaleonScripts = os.path.dirname(temp)

''' Argument parser definition '''
def ValidateARGS(args):
    if args.EvolutionaryDistancesUse == "enabled" and args.ClusterAnalysisType == "BetweenFamilies":
        raise argparse.ArgumentError(None, "This combination is not allowed: '-e enabled' and '-F BetweenFamilies' values.")

    elif args.ProteinDirectory != None and args.ClusterAnalysisType == "BetweenFamilies":
        raise argparse.ArgumentError(None, "This combination is not allowed: '-p dir_name' and '-F BetweenFamilies' values.")
    
    elif args.ProteinDirectory != None and args.EvolutionaryDistancesUse == "disabled":
        raise argparse.ArgumentError(None, "Enable the use of Evolutionary distances: '-e enabled'.")

    elif args.EvolutionaryDistancesUse == "enabled" and args.ProteinDirectory == None:
        raise argparse.ArgumentError(None, "The use of 'Evolutionary distances' is enabled. Please provide the protein directory: '-p ProtDir_name'.")

    elif os.path.exists(args.AnnotationDirectory) == False:
        raise argparse.ArgumentError(None, f"The input '--AnnotationDirectory' doesn't exist: '-a {args.AnnotationDirectory}'.")

    elif args.EvolutionaryDistancesUse == "enabled" and os.path.exists(args.ProteinDirectory) == False:
        raise argparse.ArgumentError(None, f"The input '--ProteinDirectory' doesn't exist: '-p {args.ProteinDirectory}'.")

    else:
        pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='General pipeline', formatter_class=argparse.RawTextHelpFormatter) #ArgumentDefaultsHelpFormatter)
    parser._optionals.title = "Input arguments"

    ### ~~~~ Create sub-parser ~~~~ ###

    sub_parsers = parser.add_subparsers(
        title = "Operating modes",
        dest = "mode",
        metavar = "-m",
        required = True)


    ### ~~~~ Define subparsers for each mode ~~~~ ###

    ''' gestimate - Compute several statistics for one or more gvalues '''

    parser_gEstimation = sub_parsers.add_parser("gestimate", help = "Estimate g values")


    parser_gEstimation.add_argument("-n", "--NumOfGenes",
        type = int,
        required = True,
        help = "(REQUIRED). Gene family size - number of gene copies (i.e.: 411).")

    parser_gEstimation.add_argument("-s", "--GenomeSize",
        type = float,
        required = True,
        help = "(REQUIRED). Genome size in Mb units (i.e.: 1365.69).")
    
    parser_gEstimation.add_argument("-g", "--gValue",
        nargs = "+",
        required = True,
        help = "(REQUIRED). Specify the maximum distance (in Kb units) between two copies of a given family to consider that they are clustered (i.e.: 100 -> 100 Kb). A single value or a list of comma separated values can be provided as input.")    
    
    parser_gEstimation.add_argument("-scripts", "--ScriptDirectory",
        default = "Scripts",
        type = str,
        help = "Directory with pipeline scripts. (Default: Scripts)")

    parser_gEstimation.add_argument("-log", "--LogDir",
        default = "Logs_gestimate_mode",
        type = str,
        help = "Directory with logs from the main scripts. (Default: Logs_gestimate_mode)")

    parser_gEstimation.add_argument("-outdir", "--OutDir",
        default = "g_estimation_Results_Directory",
        type = str,
        help = "Results directory for gestimate mode. (Default: g_estimation_Results_Directory)")

    ''' default - Use Physical distances and Evolutionary distance '''

    parser_default = sub_parsers.add_parser("clusterfinder", help = "Default mode (Clustering based on Physical and Evolutionary distance)")

    # Add arguments
    parser_default.add_argument("-a", "--AnnotationDirectory",
        type = str,
        required = True,
        help = "(REQUIRED). Input directory with annotation files")
    
    
    parser_default.add_argument("-p", "--ProteinDirectory",
        type = str,
        help = "Input directory with protein files (*fasta). (REQUIRED if '--EvolutionaryDistancesUse' is enabled)")

    parser_default.add_argument("-pm", "--ProteinMultipleAlignment",
        type = str,
        default = "False",
        choices = ["True", "False"],
        help = "Input directory with protein MSA files (*aln). (REQUIRED if '--EvolutionaryDistancesUse' is enabled)")

    # Optional arguments (configured with some default settings)

    parser_default.add_argument("-F", "--ClusterAnalysisType",
        default="WithinFamilies",
        type = str,
        choices = ["WithinFamilies", "BetweenFamilies"],
        help = "'Within Families' - each family is analyzed separately. 'Between families' - cluster search will be carried between two input families")

    parser_default.add_argument("-r", "--ExpansionRounds",
        default="inf",
        type = str,
        help = "Number of iterative rounds to expand the clusters. (Default: 'inf' or an integer (i.e.: 1)). Note that the 'cluster expansion' usually ends after 2 rounds")


    parser_default.add_argument("-g", "--gValue",
        default = 100,
        nargs = "+",
        help = "Specify the maximum distance between two copies of a given family to consider that they are clustered in a unit of Kb (i.e.: 100). A list of g values is also accepted: (i.e.: 0.5,10,50,100,120.5)")


    parser_default.add_argument("-scripts", "--ScriptDirectory",
        default = "Scripts",
        type = str,
        help = "Directory with pipeline scripts")

    ##############################################################################

    parser_default.add_argument("-e", "--EvolutionaryDistancesUse",
        default = "disabled",
        type = str,
        choices = ["enabled", "disabled"],
        help = "Enable or Disable the use of Evolutionary distances information in Heatmaps")


    parser_default.add_argument("-emx_pos", "--EvoMatrixPositionOnHeatmap",
        default = "Upper",
        type = str,
        choices = ["Lower", "Upper"],
        help = "Position of the Evo. distance semi-matrix on a Heatmap")

    parser_default.add_argument("-emx_sep", "--EvoMatrixFieldSeparator",
        default = "tab",
        type = str,
        choices = ["tabular", "comma", "space"],
        help = "Admitted data field separators for Evo. dist. matrices: 'tabular (\\t)', 'comma (,)', 'single space ( )'")

    ##############################################################################
    parser_default.add_argument("-c", "--ColorScaleOption",
        default = "one",
        type = str,
        choices = ["one", "two"],
        help = "Choose 'one' or 'two' color maps to be used for visualization of Physical and Evolutionary distances on a merged distance heatmap'")

    parser_default.add_argument("-f", "--SquareFrameColor",
        default = "black",
        type = str,
        help = "Choose the color of square frame that will highlight clustered genes")

    parser_default.add_argument("-f2", "--SquareFrameColor2",
        default = "blue",
        type = str,
        help = "Choose the color of square frame that will highlight clustered genes (enabled in the 'joint analysis' of two gene families to highlight clusters of the second family)")

    parser_default.add_argument("-f3", "--SquareFrameColor3",
        default = "red",
        type = str,
        help = "Choose the color of square frame that will highlight clustered genes between families")

    parser_default.add_argument("-cmap_1", "--ColorMap_1",
        default = "default",
        type = str,
        help = "Color palette for Physical Distances. Default palette: 'mako' ")
   
    parser_default.add_argument("-cmap_2", "--ColorMap_2",
        default = "default",
        type = str,
        help = "Color palette for Evolutionary Distances. Default palette: 'viridis'")

    parser_default.add_argument("-l", "--LegendScaleUnits",
        default = "default",
        choices = ["default", "auto"],
        type = str,
        help = "The heatmap's scale units format is set by default, but it can be adjusted automatically")

    parser_default.add_argument("-edec", "--EvoScale_round_decimals",
        default = "2",
        type = str,
        help = "Evolutionary scale distance round decimals on merged distance heatmap")
    
    parser_default.add_argument("-pdec", "--PhysScale_round_decimals",
        default = "2",
        type = str,
        help = "Physical scale distance round decimals on merged distance heatmap")

    ##############################################################################
    parser_default.add_argument("-b", "--Binaries",
        default = "bin",
        type = str,
        help = "PATH to bin folder containing the required programs")

    parser_default.add_argument("-bedtools_path", "--Bedtools",
        default = "bin/bedtools",
        type = str,
        help = "PATH to the bedtools binary file")

    parser_default.add_argument("-t", "--EvoTreeSoft",
        default = "FastTree",
        type = str,
        choices = ["FastTree", "iqtree"],
        help = "Software for the computation of Evolutionary Distances. (Default: FastTree)")

    parser_default.add_argument("-m", "--MSAalgorithm",
        default = "mafft-linsi",
        type = str,
        choices = ["mafft-linsi", "mafft"],
        help = "Program to conduct the multiple sequence alignment. Specify 'mafft-linsi' or 'mafft', mafft-linsi is more accurate but slow for hundreds of protein sequences (Default=mafft-linsi)")

    ##############################################################################


    ##############################################################################
    parser_default.add_argument("-o", "--PhysicalMatricesOutDir",
        default = "PhysicalDist_Matrices",
        type = str,
        help = "Output directory for Physical distance matrices. (Default: PhysicalDist_Matrices)")


    parser_default.add_argument("-log", "--LogDir",
        default = "Logs_clusterfinder_mode",
        type = str,
        help = "Directory with logs from the main scripts. (Default: Logs_clusterfinder_mode)")

    parser_default.add_argument("-outdir", "--OutDir",
        default = "clusterfinder_Results_Directory",
        type = str,
        help = "Results directory for clusterfinder mode. (Default: clusterfinder_Results_Directory)")

    ##############################################################################

    # Parsing arguments
    args = parser.parse_args()


    # Get configvalue
    config = vars(args)
    print(config)

    # Validate function
    if config["mode"] == "gestimate":
        pass
    else:
        ValidateARGS(args)

''' Functions '''

def cr_dir(i_dir): 
    # Fx to create/overwrite a directory
    
    if os.path.isdir(i_dir):
        # print(f"Removing old directory and creating a new one: '{i_dir}'")
        shutil.rmtree(i_dir)
        os.mkdir(i_dir)
    
    else:
        # print(f"Creating a new directory: '{i_dir}'")
        os.mkdir(i_dir)

def print_header(i_opt, i_text=None):
    if i_opt == "main":
        print("\n")
        print("-"*120)
        print(i_text)
        print("-"*120)
    if i_opt == "sub":
        print("~"*120)

''' Execution '''
# Create a Log directory
Log_dirname = config["LogDir"]



''' gestimate mode '''

if config["mode"] == "gestimate":
    print(f"# SELECTED MODE: {config['mode']}")
    
    # Print configuration
    for opt, opt_value in config.items():
        print(f"    {opt}   {opt_value}")

    # Create out log directory
    # print("\n# Creating LogDir")
    cr_dir(Log_dirname)    
    

    # Set the parameteres for the Rscript
    gene_number_r, genome_size_r, g_values_data = str(config["NumOfGenes"]), str(config["GenomeSize"]), config["gValue"]

    # Pre-processing of gvalue(s) for Rscript
    if len(g_values_data) == 1:
        g_values_data = str(g_values_data[0])

    else:
        temp_list = list(map(str, g_values_data))
        g_values_data = ",".join(temp_list)

    # Create out results directory
    Results_dir = config["OutDir"]
    cr_dir(Results_dir)

    # Script location
    Script_dir = config["ScriptDirectory"]
    Script_dir = f"{PATHtoGaleonScripts}/{Script_dir}"
    gestimation_rscript = f"{Script_dir}/Rscript_estimate_g.R"

    # Log files
    g_stdout_file = f"{Log_dirname}/gestimation.out"
    g_stderr_file = f"{Log_dirname}/gestimation.err"


    # Start...
    print(f"\n# Running '{gestimation_rscript}' script")

    with open(g_stdout_file, "w") as g_out, open(g_stderr_file, "w") as g_err:
        print("## cmd: ", " ".join(["Rscript", gestimation_rscript, gene_number_r, genome_size_r, g_values_data, Results_dir]))
    
        temp = subprocess.run(["Rscript", gestimation_rscript, gene_number_r, genome_size_r, g_values_data, Results_dir], text=True, capture_output=True)
        g_out.write(temp.stdout)
        g_err.write(temp.stderr)

    print("\nFinish.")
    print("Check the results:")
    print("- Estimation table: g_estimation.table.txt")
    print("- Log files: ", g_stdout_file, g_stderr_file)

''' default mode '''

if config["mode"] == "clusterfinder":
    print(f"# SELECTED MODE: {config['mode']}")

    # Print configuration
    for opt, opt_value in config.items():
        print(f"    {opt}   {opt_value}")

    # Create out log directory
    print(f"\n# Creating LogDir: {Log_dirname}")
    cr_dir(Log_dirname)  

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ### S1. Pipeline
    ### - Multiple Sequence Alignment
    ### - Tree inference
    ### - Evolutionary Distance Estimation ###

    # Scripts and Binary directories

    Script_dir = config["ScriptDirectory"]
    Script_dir = f"{PATHtoGaleonScripts}/{Script_dir}"
    
    Binaries_dir = config["Binaries"]
    Binaries_dir = f"{PATHtoGaleonScripts}/{Binaries_dir}"

    Alignment_algorithm = config["MSAalgorithm"]
    Bedtools_dir = config["Bedtools"]

    current_dir = os.getcwd()

    # Log files
    Log_dirname = config["LogDir"]

    Results_dirname = config["OutDir"]
    cr_dir(Results_dirname)

    # Main scripts
    EvoDist_exec_script = f"{Script_dir}/S1_TreeInference_and_EvoDistEstimation.sh"
    EvoDist_master_script = f"{Script_dir}/run_evodistance.sh"
    EvoDist_exec_script_pm = f"{Script_dir}/S1_TreeInference_and_EvoDistEstimation_pm.sh" # pm: when using precomuted protein MSA
    EvoDist_master_script_pm = f"{Script_dir}/run_evodistance_fromaln.sh" # pm: when using precomuted protein MSA

    # Number of analyzed families
    FamilyNum = config["ClusterAnalysisType"]
    ExpRoundParam = config["ExpansionRounds"]

    # Annotation directory
    Annot_dir = config["AnnotationDirectory"]
    Annot_dir = f"{current_dir}/{Annot_dir}"

    # Output Physical matrix directory
    PMatrix_dir = config["PhysicalMatricesOutDir"]
    PMatrix_dir = f"{Results_dirname}/{PMatrix_dir}"

    # Use of Evolutionary distances
    EvoDistUse_opt = config["EvolutionaryDistancesUse"]
    EMatrix_position = config["EvoMatrixPositionOnHeatmap"]
    EMatrix_separator = config["EvoMatrixFieldSeparator"]

    # Heatmap
    HColorScaleOpt = config["ColorScaleOption"]
    HSquareFrameColor = config["SquareFrameColor"]
    HSquareFrameColor2 = config["SquareFrameColor2"]
    HSquareFrameColor3 = config["SquareFrameColor3"]

    # HCmap_uniq = config["ClustColorMap"]
    HCmap1 = config["ColorMap_1"]
    if HCmap1 == "default":
        HCmap1 = "mako"
    HCmap2 = config["ColorMap_2"]
    if HCmap2 == "default":
        HCmap2 = "jet_r"

    # Legend scale units
    LegendScaleUnit = config["LegendScaleUnits"]
    Scale1_decimals = config["EvoScale_round_decimals"]
    Scale2_decimals = config["PhysScale_round_decimals"]

    # Protein directory
    Input_prot_dir = config["ProteinDirectory"]
    Input_prot_msa = config["ProteinMultipleAlignment"]
    EMatrix_location = Input_prot_dir

    # Set the parameters for the Evo. Dist. computation
    Tree_Soft = config["EvoTreeSoft"]
    g_values_data = config["gValue"]

    # Start...
    print("\n")


    if EvoDistUse_opt == "enabled":
        if os.path.exists(Input_prot_dir) == False:
            msg = f"The input Protein directory '{Input_prot_dir}' doesn't exist!!!"
            raise ValueError(msg)

        print_header("main", "# Starting S1. Pipeline: Multiple Sequence Alignment, Tree Inference and Evolutionary Distance Estimation #")

        # S1. out log files
        s1_stdout_file = f"{Log_dirname}/S1.out"
        s1_stderr_file = f"{Log_dirname}/S1.err"

        with open(s1_stdout_file, "w") as s1out, open(s1_stderr_file, "w") as s1err:
            
            if Input_prot_msa == "False":
                print("## cmd: ", " ".join(["bash", EvoDist_exec_script, Input_prot_dir, Tree_Soft, EvoDist_master_script, Binaries_dir]))
                print("\n")

            elif Input_prot_msa == "True":
                print("## cmd: ", " ".join(["bash", EvoDist_exec_script_pm, Input_prot_dir, Tree_Soft, EvoDist_master_script_pm, Binaries_dir]))
                print("\n")

            
            if os.path.exists(f"{Binaries_dir}/{Tree_Soft}"):
                if os.access(f"{Binaries_dir}/{Tree_Soft}", os.X_OK):
                    try:
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                        if Input_prot_msa == "False":
                            temp = subprocess.run(["bash", EvoDist_exec_script, Input_prot_dir, Tree_Soft, EvoDist_master_script, Binaries_dir], text=True, capture_output=True, check=True)
                            print(temp.stdout)

                            s1out.write(temp.stdout)
                            s1err.write(temp.stderr)

                        elif Input_prot_msa == "True":
                            temp = subprocess.run(["bash", EvoDist_exec_script_pm, Input_prot_dir, Tree_Soft, EvoDist_master_script_pm, Binaries_dir], text=True, capture_output=True, check=True)
                            print(temp.stdout)

                            s1out.write(temp.stdout)
                            s1err.write(temp.stderr)

                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    except subprocess.CalledProcessError as e:
                        # raise ValueError(e.stdout)
                        raise ValueError(e.stderr)
                else:
                    emsg = f"Error! Tree Inference Software: '{Binaries_dir}/{Tree_Soft}' is not executable. Make it executable with 'chmod +x'"
                    raise ValueError(emsg)

            else:
                emsg = f"Error! Tree Inference Software: '{Binaries_dir}/{Tree_Soft}' doesn't exist"
                raise ValueError(emsg)

        
        print("\nFinish.")
        print("Check the results:")
        print(f"- Evo. Distance matrices (*distancematrix.tsv) in '{Input_prot_dir}' directory")
        print("- Log files: ", s1_stdout_file, s1_stderr_file)

    elif EvoDistUse_opt == "disabled":
        print_header("main", "# Skipping S1... #")

    else:
        raise ValueError()

    ### S2. Pipeline
    ### - Cluster definition based on Physical distances.
    ### - Inclusion of Evolutionary distance information
    ### - Mann-Whitney analysis

    if FamilyNum == "WithinFamilies":
        # Scripts
        S21_script = f"{Script_dir}/S21_ProcessInputData.py"
        S22_script = f"{Script_dir}/S22_ProcessGFF.py"
        S23_script = f"{Script_dir}/S23_Searchclusters.py"
        S24_script = f"{Script_dir}/S24_Heatmaps.py"
        S25_script = f"{Script_dir}/S25_ScatterPlots.py"

        # Start...
        print_header("main", "# Starting S2. Pipeline: Cluster Definition #")

        # S21
        print("\n")
        print_header("sub")
        print("## S2.1 Annotation Files Processing #\n")

        s21_file_log = f"{Log_dirname}/S21.log"
            
        with open(s21_file_log, "w") as s21out:
            if "bedtools" in os.listdir(Binaries_dir):
                Bedtools_dir = f"{Binaries_dir}/bedtools"
            else:
                pass

            print("## cmd: ", " ".join(["python", S21_script, Annot_dir, Bedtools_dir]))

        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            try:
                temp = subprocess.run(["python", S21_script, Annot_dir, Bedtools_dir], text=True, capture_output=True, check=True)
                print(temp.stdout)

                s21out.write(temp.stdout)
                s21out.write(temp.stderr)

            except subprocess.CalledProcessError as e:
                emsg = f"Script: {S21_script} failed. Please check that: \
                \n 1) The bedtools was correctly provided \
                \n - PATH to the bedtools binary file: '{Bedtools_dir}' \
                \n 2) The input files are in ONLY one of the following formats: (gff3, bed1, bed2) and end like that: (*_fam.gff3, *_fam.bed1, *_fam.bed2) \
                \n - PATH to the input files: '{Annot_dir}'"

                raise ValueError(emsg)


        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # S22
        print_header("sub")
        print("## S2.2 Compute Pairwise Physical Distances #\n")

        s22_file_log = f"{Log_dirname}/S22.log"

        with open(s22_file_log, "w") as s22out:

            print("## cmd: ", " ".join(["python", S22_script, Annot_dir, PMatrix_dir]))

            # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            temp = subprocess.run(["python", S22_script, Annot_dir, PMatrix_dir], text=True, capture_output=True, check=True)

            print(temp.stdout)

            s22out.write(temp.stdout)
            s22out.write(temp.stderr)

            # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # S23
        print_header("sub")
        print("## S2.3 Clusters search #")

        if isinstance(g_values_data, int):
            g_values_list = [int(g_values_data)]

        elif isinstance(g_values_data, float):
            g_values_list = [int(g_values_data)]

        elif isinstance(g_values_data, list):
            temp = g_values_data[0].split(",")
            g_values_list = [float(tt) for tt in temp]

        else:
            raise ValueError("Unknown input format for 'g_values_data'")
        
        # Run scripts for each input g value
        for g_val in g_values_list:
            g_val = float(g_val)

            print(f"\n=> Input g = {g_val}\n")
            
            log_file_3 = f"{Log_dirname}/S23.g={g_val}kb.log"
            
            with open(log_file_3, "w") as l3:

                print("## cmd: ", " ".join(["python", S23_script, str(g_val), ExpRoundParam, PMatrix_dir, FamilyNum]))

            # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                subprocess.run(["python", S23_script, str(g_val), ExpRoundParam, PMatrix_dir, FamilyNum], stdout=l3) # Eliminar el Scale value i deixar 1e5 en el propi script
            # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            # S24
            print_header("sub")
            print("## Visualization #")
            print_header("sub")
            print("## S2.4 Heatmaps #")
                    
            FamDir_list = [f for f in os.listdir(PMatrix_dir) if f.endswith("_matrices")]
            print(FamDir_list)
            log_file_4 = f"{Log_dirname}/S24.g={g_val}kb.log"
            log_file_stats = f"{Log_dirname}/S24.stats.g={g_val}kb.log"

            # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            with open(log_file_4, "w") as l4:
                for FamName in FamDir_list:
                    # print(FamName)

                    
                    if EvoDistUse_opt == "enabled":

                        print("## cmd: "," ".join(["python", 
                                        
                                        S24_script, 
                                        EvoDistUse_opt,
                                        str(g_val), 
                                        FamName, PMatrix_dir, 
                                        EMatrix_location, EMatrix_position, EMatrix_separator, LegendScaleUnit, Scale1_decimals, Scale2_decimals,
                                        HColorScaleOpt, HSquareFrameColor, HCmap1, HCmap2]))

                        Dparam = {"S24 script name: ": S24_script,
                                            "g value: ": str(g_val),
                                            "Family Name: ": FamName,
                                            "Physical Matrix directory: ": PMatrix_dir,
                                            "Evo. Distances Usage: ": EvoDistUse_opt,
                                            "Evo. Matrix location: ": EMatrix_location,
                                            "Evo. Matrix position: ": EMatrix_position,
                                            "Evo. Matrix separator: ": EMatrix_separator,
                                            "Legend scale units: ": LegendScaleUnit,
                                            "Evo. dist. legend scale decimals: ": Scale1_decimals,
                                            "Phys. dist. legend scale decimals: ": Scale2_decimals,
                                            "Color scale option: ": HColorScaleOpt, 
                                            "Square frame color: "  : HSquareFrameColor,
                                            "Color map 1: " : HCmap1,
                                            "Color map 2: " : HCmap2}
                        
                        # Print input parameteres to S24 script                                        
                        for Dk, Dvalue in Dparam.items():
                            print("- ", Dk, Dvalue)

                        # Create heatmaps
                        subprocess.run(["python", 
                                        
                                        S24_script, 
                                        EvoDistUse_opt,
                                        str(g_val), 
                                        FamName, PMatrix_dir, 
                                        EMatrix_location, EMatrix_position, EMatrix_separator, LegendScaleUnit, Scale1_decimals, Scale2_decimals,
                                        HColorScaleOpt, HSquareFrameColor, HCmap1, HCmap2],                          
                                        stdout=l4)

                        

                        # Create scatterplots
                        print_header("sub")
                        print("## S2.5 ScatterPlots #")


                        log_file_5 = f"{Log_dirname}/S25.g={g_val}kb.log"
                
                        with open(log_file_5, "w") as l5:

                            print("## cmd: ", " ".join(["python", 
                                        
                                        S25_script, 
                                        str(g_val), 
                                        FamName, PMatrix_dir, 
                                        "MergedDistances_Dataframes"]))

                            subprocess.run(["python", 
                                            
                                        S25_script, 
                                        str(g_val), 
                                        FamName, PMatrix_dir, 
                                        "MergedDistances_Dataframes"],                                
                                        stdout=l5)

                    elif EvoDistUse_opt == "disabled":

                        print("## cmd: ", " ".join(["python", 
                                                                
                                        S24_script, 
                                        EvoDistUse_opt,
                                        str(g_val), 
                                        FamName, PMatrix_dir, FamilyNum,
                                        HColorScaleOpt, HSquareFrameColor, HSquareFrameColor2, HSquareFrameColor3, HCmap1, HCmap2]))

                        Dparam = {"S24 script name: ": S24_script,
                                        "Evo. Distances Usage: ": EvoDistUse_opt,
                                        "g value: ": str(g_val),
                                        "Family Name: ": FamName,
                                        "Physical Matrix directory: ": PMatrix_dir,
                                        "Color scale option: ": HColorScaleOpt, 
                                        "Square frame color 1: "  : HSquareFrameColor,
                                        "Color map 1: " : HCmap1}

                        # Print input parameteres to S24 script
                        for Dk, Dvalue in Dparam.items():
                            print("- ", Dk, Dvalue)

                        subprocess.run(["python", 
                                        S24_script, 
                                        EvoDistUse_opt,
                                        str(g_val), 
                                        FamName, PMatrix_dir, FamilyNum, 
                                        HColorScaleOpt, HSquareFrameColor, HSquareFrameColor2, HSquareFrameColor3, HCmap1, HCmap2],
                                        stdout=l4)
                    
            # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    elif FamilyNum == "BetweenFamilies":
        
        # Scripts
        S21_script = f"{Script_dir}/S21_ProcessInputData_BetweenFamilies.py"
        S22_script = f"{Script_dir}/S22_ProcessGFF_BetweenFamilies.py"
        S23_script = f"{Script_dir}/S23_Searchclusters.py"
        S24_script = f"{Script_dir}/S24_Heatmaps.py"

        # Start...
        print_header("main", "# Starting S2. Pipeline: Cluster Definition for Two Families #")

        # S21
        print("\n")
        print_header("sub")
        print("## S2.1 Annotation Files Processing #\n")

        s21_file_log = f"{Log_dirname}/S21.TwoFam.log"
            
        with open(s21_file_log, "w") as s21out:
            if "bedtools" in os.listdir(Binaries_dir):
                Bedtools_dir = f"{Binaries_dir}/bedtools"
            else:
                pass

            print("## cmd: ", " ".join(["python", S21_script, Annot_dir, Bedtools_dir]))

            try:
                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                temp = subprocess.run(["python", S21_script, Annot_dir, Bedtools_dir], text=True, capture_output=True, check=True)
                print(temp.stderr)
                s21out.write(temp.stdout)
                s21out.write(temp.stderr)

                #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            except subprocess.CalledProcessError as e:
                # raise ValueError(e.stdout)
                raise ValueError(e.stderr)

        # S22
        print_header("sub")
        print("## S2.2 Compute Pairwise Physical Distances #\n")

        s22_file_log = f"{Log_dirname}/S22.TwoFam.log"

        with open(s22_file_log, "w") as s22out:

            print("## cmd: ", " ".join(["python", S22_script, Annot_dir, PMatrix_dir]))

            # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            temp = subprocess.run(["python", S22_script, Annot_dir, PMatrix_dir], text=True, capture_output=True, check=True)

            s22out.write(temp.stdout)
            s22out.write(temp.stderr)

            # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # S23
        print_header("sub")
        print("## S2.3 Clusters search #")

        if isinstance(g_values_data, int):
            g_values_list = [int(g_values_data)]

        elif isinstance(g_values_data, float):
            g_values_list = [int(g_values_data)]

        elif isinstance(g_values_data, list):
            temp = g_values_data[0].split(",")
            g_values_list = [float(tt) for tt in temp]

        else:
            raise ValueError("Unknown input format for 'g_values_data'")
        
        # Run scripts for each input g value
        for g_val in g_values_list:
            g_val = float(g_val)

            print(f"\n=> Input g = {g_val}\n")
            
            log_file_3 = f"{Log_dirname}/S23.g={g_val}kb.TwoFam.log"
            
            with open(log_file_3, "w") as l3:

                print("## cmd: ", " ".join(["python", S23_script, str(g_val), ExpRoundParam, PMatrix_dir, FamilyNum]))
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                subprocess.run(["python", S23_script, str(g_val), ExpRoundParam, PMatrix_dir, FamilyNum], stdout=l3) # Eliminar el Scale value i deixar 1e5 en el propi script
            #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            # S24
            print_header("sub")
            print("## Visualization #")
            print_header("sub")
            print("## S2.4 Heatmaps #")
                    
            FamDir_list = [f for f in os.listdir(PMatrix_dir) if f.endswith("_matrices")]
            
            log_file_4 = f"{Log_dirname}/S24.g={g_val}kb.TwoFam.log"

            # #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            with open(log_file_4, "w") as l4:
                for FamName in FamDir_list:
                    print(FamName)

                    if EvoDistUse_opt == "disabled":

                        print("## cmd: ", " ".join(["python", 
                                        
                                        S24_script, 
                                        EvoDistUse_opt,
                                        str(g_val), 
                                        FamName, PMatrix_dir,
                                        FamilyNum, HColorScaleOpt, 
                                        HSquareFrameColor, HSquareFrameColor2, HSquareFrameColor3, HCmap1, HCmap2]))

                        Dparam = {"S24 script name: ": S24_script,
                                        "Evo. Distances Usage: ": EvoDistUse_opt,
                                        "g value: ": str(g_val),
                                        "Family Name: ": FamName,
                                        "Physical Matrix directory: ": PMatrix_dir,
                                        "Color scale option: ": HColorScaleOpt, 
                                        "Square frame color 1: "  : HSquareFrameColor,
                                        "Square frame color 2: "  : HSquareFrameColor2,
                                        "Square frame color 3: "  : HSquareFrameColor3,
                                        "Color map 1: " : HCmap1,
                                        "Color map 2: " : HCmap2}
                        
                        # Print input parameteres to S24 script
                        for Dk, Dvalue in Dparam.items():
                            print("- ", Dk, Dvalue)

                        subprocess.run(["python", 
                                        
                                        S24_script, 
                                        EvoDistUse_opt,
                                        str(g_val), 
                                        FamName, PMatrix_dir, FamilyNum, HColorScaleOpt,
                                        HSquareFrameColor, HSquareFrameColor2, HSquareFrameColor3, HCmap1, HCmap2], 
                                        stdout=l4)
                    else:
                        raise ValueError()
