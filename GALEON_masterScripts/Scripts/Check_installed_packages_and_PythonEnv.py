import argparse, sys

# Argument Parser
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check script', formatter_class=argparse.RawTextHelpFormatter) #ArgumentDefaultsHelpFormatter)
    # parser._optionals.title = "Input arguments"
    ## ~~~~ Create sub-parser ~~~~ ###

    parser.add_argument( "-bedtools", "--BedtoolsPATH",
        type = str,
        default = "bin/bedtools",
        help = "PATH to bedtools")

    parser.add_argument( "-fasttree", "--FastTreePATH",
        type = str,
        default = "bin/FastTree",
        help = "PATH to FastTree")

    parser.add_argument( "-iqtree", "--IqTree2PATH",
        type = str,
        default = "bin/iqtree2",
        help = "PATH to iqtree2")

    parser.add_argument( "-mafft", "--MafftPATH",
        type = str,
        default = "mafft",
        help = "PATH to mafft")

    # Parsing arguments
    args = parser.parse_args()
    
    # Get configvalue
    config = vars(args)
    print(config)

D_SOFTNAMES = {"BedtoolsPATH" : "bedtools", 
                "FastTreePATH" : "FastTree", 
                "MafftPATH" : "mafft", 
                "IqTree2PATH" : "iqtree2"}

D_SOFT2PATH = {"bedtools" : "BedtoolsPATH", 
                "FastTree" : "FastTreePATH", 
                "mafft" : "MafftPATH", 
                "iqtree2" : "IqTree2PATH"}


# Check Python version
if sys.version_info < (3, 0):
    sys.exit("This script requires Python 3 or higher. You are using Python 2.x.")

packages_to_check = [
    'argparse',
    'ast',
    'collections',
    'copy',
    'gc',
    'itertools',
    'matplotlib',
    'numpy',
    'operator',
    'os',
    'pandas',
    're',
    'seaborn',
    'shutil',
    'string',
    'scipy',
    'subprocess',
    'sys',
    'time',
]

missing_packages = []

for package in packages_to_check:
    try:
        __import__(package)
    except ImportError:
        missing_packages.append(package)

if missing_packages:
    print("The following packages/programs are missing:")
    for item in missing_packages:
        print(item)
else:
    print("All the packages and binaries are installed and available.")


import os, subprocess
# Export bin directory
BinDir = "bin/"
SoftList = ["FastTree", "iqtree2", "bedtools"]

if os.path.exists(BinDir) == False:
    print(f"No {BinDir} directory detected. \
        \n- Trying to get software paths from the arguments: -bedtools, -fasttree, -iqtree")

    for kSOFT, vPATH in config.items():
        if os.path.exists(vPATH) == False:
            emsg = f"{kSOFT}. PATH not found: '{vPATH}'"
            raise ValueError(emsg)
        else:
            softName, softPath = kSOFT, os.path.abspath(vPATH)
            print(f"{softName} : {softPath}")

            softPath = f"{softPath}/{D_SOFTNAMES[softName]}"

            # Check if the binary file exists
            if os.path.exists(softPath) :
                msg = f"- Using '{D_SOFTNAMES[softName]}' from {softPath}"
                print(msg)
                
                # Check if the binary file is executable
                if os.access(softPath, os.X_OK):
                    # Get the current PATH and append the binary's directory to it
                    os.environ['PATH'] = f"{os.path.dirname(softPath)}:{os.environ['PATH']}"
                else:
                    emsg = f"This binary file '{softPath}' is not executable"
                    raise ValueError(emsg)
            else:
                emsg = f"Something is wrong, this binary file '{softPath}' doesn't exist"
                raise ValueError(emsg)
else:
    print(f"{BinDir} directory detected.")

    for soft in SoftList:
        if soft not in os.listdir(BinDir):
            softPath = config[D_SOFT2PATH[soft]]
            if os.path.exists(softPath):
                msg = f"- Using '{soft}' from {softPath}"
                # Check if the binary file is executable
                if os.access(softPath, os.X_OK):
                    print(msg)
                    # Get the current PATH and append the binary's directory to it
                    os.environ['PATH'] = f"{os.path.dirname(softPath)}:{os.environ['PATH']}"
                else:
                    emsg = f"This binary file '{softPath}' is not executable"
                    raise ValueError(emsg)
            else:
                emsg = f"Mandatory software '{softPath}' does not exist"
                raise ValueError(emsg)
        else:
            softPath = os.path.abspath(BinDir)
            msg = f"- Using '{soft}' from {softPath}"
            print(msg)

            # Check if the binary file is executable
            softPath = f"{softPath}/{soft}"
            if os.access(softPath, os.X_OK):
                # Get the current PATH and append the binary's directory to it
                os.environ['PATH'] = f"{os.path.dirname(softPath)}:{os.environ['PATH']}"
            else:
                emsg = f"This binary file '{softPath}' is not executable"
                raise ValueError(emsg)            

# Check mafft
temp = subprocess.run("which mafft", shell=True, capture_output=True, text=True).stdout.strip()
if temp == "":
    raise ValueError("'mafft' not found, check whether 'mafft' is installed")
else:
    msg = f"- Using 'mafft' from {temp}"
    print(msg)


# Check newick_utils
temp = subprocess.run("which nw_distance", shell=True, capture_output=True, text=True).stdout.strip()
if temp == "":
    raise ValueError("'nw_distance' not found, check whether 'newick_utils' is installed")
else:
    msg = f"- Using 'nw_distance' from {temp}"
    print(msg)
    
# Check R
temp = subprocess.run("which Rscript", shell=True, capture_output=True, text=True).stdout.strip()
if temp == "":
    raise ValueError("'Rscript' not found, check whether 'Rscript' is installed")
else:
    msg = f"- Using 'Rscript' from {temp}"
    print(msg)


# Check R
temp = subprocess.run("which R", shell=True, capture_output=True, text=True).stdout.strip()
if temp == "":
    raise ValueError("'R' not found, check whether 'R' is installed")
else:
    msg = f"- Using 'R' from {temp}"
    print(msg)


# Check pandoc
temp = subprocess.run("which pandoc", shell=True, capture_output=True, text=True).stdout.strip()
if temp == "":
    raise ValueError("'pandoc' not found")
else:
    msg = f"- Using 'pandoc' from {temp}"
    print(msg)
