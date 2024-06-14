#%%
import subprocess
import os


# enter the path to gcc and g++ compilers
gcc = "D:/msys2/mingw64/bin/gcc.exe"
gpp = "D:/msys2/mingw64/bin/g++.exe"

# for windows
powershell = "C:/Windows/System32/WindowsPowerShell/v1.0/powershell.exe"
# for Linux
# powershell = "/bin/bash"

# depends on the compiler toolchain 
generator = "MinGW Makefiles"
# generator = "Unix Makefiles"
# generator = "Ninja"
# generator = "MSYS Makefiles"


# check path 
if not os.path.exists(gcc):
    print("gcc path not found!")
    exit(1)
if not os.path.exists(gpp):
    print("g++ path not found!")
    exit(1)
if not os.path.exists(powershell):
    print("powershell/bash path not found!")
    
    
#%%
if os.path.isdir('./build/'):
    # clear previous cache 
    clear_cache_cmd = ["rm", "-r", "./build/*"]
    out = subprocess.run(clear_cache_cmd, shell=True, check=False, text=True, capture_output=True, executable=powershell)

    if out.returncode != 0:
        print(out.stdout)
        print(out.stderr)
        exit(1)
else:
    os.mkdir('./build/')
    
# config cmake
cmake_config_cmd = ["cmake", "-DCMAKE_C_COMPILER:FILEPATH="+"'{}'".format(gcc), "-DCMAKE_CXX_COMPILER:FILEPATH="+"'{}'".format(gpp), "-S", ".", "-B", "build", "-G", "'{}'".format(generator)]

out = subprocess.run(cmake_config_cmd, shell=True, check=False, text=True, capture_output=True, executable=powershell)

if out.returncode != 0:
    print(out.stdout)
    print(out.stderr)
    exit(1)

print('CMake configuration successful!')


# build 
build_command = ["cmake", "--build", "build", "--config", "Release", "--target", "all"]

out = subprocess.run(build_command, shell=True, check=False, text=True, capture_output=True, executable=powershell)

if out.returncode != 0:
    print(out.stdout)
    print(out.stderr)
    exit(1)
    
print('Build successful!')


# install 
install_command = ["cmake", "--install", "build", "--config", "Release"]

out = subprocess.run(install_command, shell=True, check=False, text=True, capture_output=True, executable=powershell)

if out.returncode != 0:
    print(out.stdout)
    print(out.stderr)
    exit(1)
    
print('Install successful!')


#%%
# run the executable
# run_command = ["./bin/OpenLPT.exe"]

# out = subprocess.run(run_command, check=True, text=True, capture_output=True)

# print(out.stderr)
# print(out.stdout)
# %%
