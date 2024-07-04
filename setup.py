#%%
import os
import subprocess
import sys
import platform
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DPYOPENLPT=True']
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
            
        out = subprocess.run(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env, text=True, capture_output=True)
        if out.returncode != 0:
            print(out.args)
            print(out.stdout)
            print(out.stderr)
            exit(1)
            
        out = subprocess.run(['cmake', '--build', '.'] + build_args, cwd=self.build_temp, text=True, capture_output=True)
        if out.returncode != 0:
            print(out.args)
            print(out.stdout)
            print(out.stderr)
            exit(1)
        
            

setup(
    name='pyOpenLPT',
    version='0.1.0',
    author='Shijie Zhong',
    author_email='szhong12@jhu.edu',
    description='A Python interface for OpenLPT',
    long_description='',
    ext_modules=[CMakeExtension('pyOpenLPT')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    install_requires=[
        'pybind11>=2.6.0',
        'numpy>=1.16.0'
    ],
    packages=find_packages(),
    include_package_data=True,
)