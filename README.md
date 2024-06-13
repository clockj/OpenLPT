# OpenLPT 2.0

OpenLPT 2.0 is a new version of particle tracking code. It can be used for tracking dense tracers, polydisperse bubbles and many other objects. 

Look how easy it is to use:
```bash 
${code_path}/bin/OpenLPT.exe config.txt
```


## Features

- Lagrangian particle tracking for multiple objects (point-like particles, spherical particles, etc.)
- Support stereomatching with multiple cameras (at least 2)
- Include multiple test cases for users to test and understand the code
- Better structure for adding new functions


## Installation

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
After finishing the above steps, OpenLPT 2.0 can be installed by running (Windows users can use **powershell** to run the following commands):

```bash
git clone https://github.com/clockj/OpenLPT.git

cd OpenLPT
mkdir build

cmake --no-warn-unused-cli -DCMAKE_BUILD_TYPE:STRING=Release -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_C_COMPILER:FILEPATH=${path_to_gcc.exe} -DCMAKE_CXX_COMPILER:FILEPATH=${path_to_g++.exe} "-S${code_path}" "-B${code_path}/build/" -G "MinGW Makefiles"

cmake --build ./build --config Release --target all

cmake --install ./build
```

Note: Depending on the gcc toolchain, the option after **-G** can also be **"Unix Makefiles"** or **Ninja**. If you are a Linux user, you probably do not need to specify this option.

In order to run the excusable program in the command line, users sometimes needs to change the permission to run scripts.

- For Windows users, you should open **powershell** as administrator. And then enter `Set-ExecutionPolicy -ExecutionPolicy RemoteSigned`. The policy can be switched back by entering `Set-ExecutionPolicy -ExecutionPolicy Restricted`


## Contribute

- Issue Tracker: 
- Source Code: [OpenLPT 2.0](https://github.com/clockj/OpenLPT.git)


## Support

If you are having issues, please let me know.
I have a mailing list located at: szhong12@jhu.edu 


## License

The project is licensed under the MIT license.

