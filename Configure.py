import os, sys, subprocess

def add_interpreter_path(directory, interpreter_path, i_language):
    if i_language == "python":
        ending = ".py"
    elif i_language == "perl":
        ending = ".pl"

    for r, d, files in os.walk(directory):
        if len(files) != 0:
            for filename in files:
                if filename.endswith(ending):
                    filepath = os.path.join(r, filename)
                    with open(filepath, 'r+') as f:
                        lines = f.readlines()
                        if not lines[0].startswith('#!'):
                            lines.insert(0, f'#!{interpreter_path}')  # Default Python interpreter
                            f.seek(0)
                            f.writelines(lines)
                            #print(f'Added shebang line to {filename}')
                        else:
                            #print(lines[0], [interpreter_path])
                            pass


# Example usage:
directory = sys.argv[1]
python_interpreter_path = subprocess.run("which python", shell=True, capture_output=True, text=True).stdout
add_interpreter_path(directory, python_interpreter_path, "python")

perl_interpreter_path = subprocess.run("which perl", shell=True, capture_output=True, text=True).stdout
add_interpreter_path(directory, perl_interpreter_path, "perl")
