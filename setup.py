# Modified by Velat Kilic from
# Source: https://github.com/pybind/cmake_example
import os
import sys
import platform
import subprocess

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.spawn import find_executable


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = sourcedir


class CMakeBuild(build_ext):
    def run(self):
        if find_executable('cmake') is None:
            print("CMake must be installed to build this extension")
            sys.exit(-1)

        for ext in self.extensions:
            build_dir = os.path.join(os.path.dirname(__file__), 'build', 'cmake')
            if not os.path.exists(build_dir):
                os.makedirs(build_dir)
            cmake_dir = os.path.abspath(ext.sourcedir)

            extpath = self.get_ext_fullpath(ext.name)
            extfulldir = os.path.abspath(os.path.dirname(extpath))
            cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extfulldir,
                          '-DPYTHON_EXECUTABLE=' + sys.executable,
                          '-DBUILD_PIP=True']
            build_args = ['--config', 'Release']

            if platform.system() == "Windows":
                cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE=' + extfulldir]
                if sys.maxsize > 2**32:
                    cmake_args += ['-A', 'x64']
                build_args += ['--', '/m']
            else:
                build_args += ['--', '-j2']

            subprocess.check_call(['cmake', cmake_dir] + cmake_args, cwd=build_dir)
            subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=build_dir)

setup(
    name='pylisa',
    version='0.0.1',
    author='Velat Kilic',
    author_email='velatkilic@gmail.com',
    description='Lidar light scattering augmentation',
    long_description='',
    ext_modules=[CMakeExtension('pylisa')],
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    keywords=['lidar', 'augmentation', 'scattering', 'lidar augmentation'],
    classifiers = [
        "Development Status :: 2 - Pre-Alpha",
    ],
)
