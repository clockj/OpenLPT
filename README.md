# OpenLPT 2.0

OpenLPT 2.0 is a new version of particle tracking code. It can be used for tracking dense tracers, polydisperse bubbles and many other objects. 

Look how easy it is to use:
```bash 
# Use in command line
${code_path}/bin/OpenLPT.exe config.txt
```

```python
# Use in python
import pyOpenLPT as lpt

# redirect std::cout to python 
redirector = lpt.PythonStreamRedirector() 

config_file = '${path_to_config_file}'
lpt.run(config_file)
```


## Features
- User-friendly interface in python
- Lagrangian particle tracking for multiple objects (point-like particles, spherical particles, etc.)
- Support stereomatching with multiple cameras (at least 2)
- Include multiple test cases for users to test and understand the code
- Better structure for adding new functions


## Installation - Python Version

### Pre-request
[CMake](https://cmake.org/),
[Anaconda](https://www.anaconda.com/) or
[Miniconda](https://docs.anaconda.com/miniconda/)

### Install OpenLPT

Create a python environment and install dependencies
```bash
# use conda
conda create -n OpenLPT 
conda activate OpenLPT
pip install -r requirements.txt
conda deactivate
```
Download the source code from github and install it
```bash
git clone https://github.com/clockj/OpenLPT.git
cd OpenLPT
conda activate OpenLPT
pip install .
conda deactivate
```
Use the package
```python
# conda activate OpenLPT (run in the bash)
import pyOpenLPT as lpt 
redirector = lpt.PythonStreamRedirector() 

config_file = '${path_to_config_file}'
lpt.run(config_file)
```

## Installation - CPP Version
Users who want to install the pure cpp version can refer to the file [README_CPP.md](README_CPP.md)

## Contribute

- Issue Tracker: 
- Source Code: [OpenLPT 2.0](https://github.com/clockj/OpenLPT.git)


## Support

If you are having issues, please let me know.
I have a mailing list located at: szhong12@jhu.edu 


## License

The project is licensed under the MIT license.

