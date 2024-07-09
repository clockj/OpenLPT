## Installation - CPP Version

### Pre-request: 
[CMake](https://cmake.org/),
[MinGW](https://www.mingw-w64.org/) (or any other gcc/g++ compilers).
**Note that the /bin/ folder of CMake and MinGW should be added to the system environment variable.**

For **Windows** users, an example to install **MinGW** is shown in the following.

1. Install [msys2](https://www.msys2.org/).
2. Open **MSYS2 MINGW64**. This software is not the one opened defaultly
3. Install the compiler toolchain:
`
pacman -S mingw-w64-x86_64-toolchain
`
4. Add the **/bin/** folder of mingw64 (`${msys2_path}/mingw64/bin/`) to the system environment variables. 



### Install OpenLPT

#### Install by python

A python code is provided to build the project. Users just need to **specify the path to gcc, g++, powershell/bash and the option of generator in the code**. 
For Windows users, if you follow the above steps, then there is no need to change the option of generator.

After that, it can be installed by running
```bash
python build.py
```

#### Install by command line

OpenLPT 2.0 can also be installed by using command line (Windows users can use **powershell** to run the following commands):

```bash
git clone https://github.com/clockj/OpenLPT.git

cd OpenLPT
mkdir build

cmake --no-warn-unused-cli -DCMAKE_BUILD_TYPE:STRING=Release -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_C_COMPILER:FILEPATH=${path_to_gcc.exe} -DCMAKE_CXX_COMPILER:FILEPATH=${path_to_g++.exe} "-S${code_path}" "-B${code_path}/build/" -G "MinGW Makefiles"

cmake --build ./build --config Release --target all

cmake --install ./build
```

Note: Depending on the compiler toolchain, the option after **-G** can also be **"Unix Makefiles"** or **Ninja**. If you are a Linux user, you probably do not need to specify this option.

In order to run the excusable program in the command line, users sometimes needs to change the permission to run scripts.

- For Windows users, you should open **powershell** as administrator. And then enter `Set-ExecutionPolicy -ExecutionPolicy RemoteSigned`. The policy can be switched back by entering `Set-ExecutionPolicy -ExecutionPolicy Restricted`