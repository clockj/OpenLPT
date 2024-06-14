#%%
import subprocess

#%%
gcc = "D:/msys2/ucrt64/bin/gcc.exe"
gpp = "D:/msys2/ucrt64/bin/g++.exe"
powershell = "C:/Windows/System32/WindowsPowerShell/v1.0/powershell.exe"
generator = "Unix Makefiles"
# generator = "Ninja"

#%%
# clear previous cache 
clear_cache_cmd = ["rm", "-r", "./build/*"]
out = subprocess.run(clear_cache_cmd, shell=True, check=False, text=True, capture_output=True, executable=powershell)

print(out.stdout)
print(out.stderr)

cmake_config_cmd = ["cmake", "-DCMAKE_C_COMPILER:FILEPATH="+gcc, "-DCMAKE_CXX_COMPILER:FILEPATH="+gpp, "-S", ".", "-B", "build", "-G", "'{}'".format(generator)]

out = subprocess.run(cmake_config_cmd, shell=True, check=False, text=True, capture_output=True, executable=powershell)

print(out.stdout)
print(out.stderr)


build_command = ["cmake", "--build", "build", "--config", "Release", "--target", "all"]

out = subprocess.run(build_command, shell=True, check=False, text=True, capture_output=True, executable=powershell)

print(out.stdout)
print(out.stderr)

#%%
run_command = ["./build/OpenLPT.exe"]

out = subprocess.run(run_command, check=True, text=True, capture_output=True)

print(out.stderr)
print(out.stdout)
# %%
