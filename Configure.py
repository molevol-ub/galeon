import os, sys, subprocess

def add_interpreter_path(directory, interpreter_path):
    for r, d, files in os.walk(directory):
        if len(files) != 0:
            for filename in files:
                if filename.endswith('.py'):
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
interpreter_path = subprocess.run("which python", shell=True, capture_output=True, text=True).stdout
add_interpreter_path(directory, interpreter_path)
