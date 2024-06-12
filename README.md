OpenLPT 2.0
========

OpenLPT 2.0 is a new version of particle tracking code. It can be used for tracking dense tracers, polydisperse bubbles and many other objects. 

Look how easy it is to use:
```bash
${code_path}/bin/OpenLPT.exe config.txt
```


Features
--------

- pending


Installation
------------

Pre-request: 
[CMake](https://cmake.org/),
[MinGW](https://www.mingw-w64.org/) (or any other gcc/g++ compilers)

After finishing the above steps, OpenLPT 2.0 can be installed by running:

```bash
git clone https://github.com/clockj/OpenLPT.git
cd OpenLPT

mkdir build

cmake --no-warn-unused-cli -DCMAKE_BUILD_TYPE:STRING=Release -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=TRUE -DCMAKE_C_COMPILER:FILEPATH=${path_to_gcc.exe} -DCMAKE_CXX_COMPILER:FILEPATH=${path_to_g++.exe} "-S${code_path}" "-B${code_path}/build/" -G "MinGW Makefiles"
```

Note: Depending on the gcc toolchain, the other options after **-G** can be **"Unix Makefiles"** or **ninja**.


Contribute
----------

- Issue Tracker: 
- Source Code: [OpenLPT 2.0](https://github.com/clockj/OpenLPT.git)


Support
-------

If you are having issues, please let us know.
We have a mailing list located at: 


License
-------

The project is licensed under the MIT license.

