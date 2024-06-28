from setuptools import setup, Extension
from pybind11.setup_helpers import Pybind11Extension, build_ext
import glob

# require pybind11, numpy

src = ["src/bindings.cpp"] + glob.glob('src/pybind_OpenLPT/*.cpp')
# src = ['src/bindings.cpp', 'src/pybind_OpenLPT/pyMatrix.cpp']
dir = glob.glob('src/*/', recursive=True) + glob.glob('inc/*/', recursive=True)

ext_modules = [
    Pybind11Extension(
        "pyOpenLPT",
        src,
        include_dirs=dir,
    ),
]

setup(
    name="pyOpenLPT",
    version="0.1.0",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
)