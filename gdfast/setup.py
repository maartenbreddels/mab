# -*- coding: utf-8 -*-
import os
from distutils.sysconfig import get_python_inc, get_python_lib
import numpy
import platform


from distutils.core import setup, Extension
numdir = os.path.dirname(numpy.__file__)

import sys 
sys.setrecursionlimit(10000)

APP = ["bin/vaex"]
DATA_FILES = []

import glob

DATA_FILES #.append(["", glob.glob("doc/*/*")])
#OPTIONS = {'argv_emulation': False, 'excludes':[], 'resources':['python/gavi/icons'], 'matplotlib_backends':'-'}



include_dirs = []
library_dirs = []
#library_dirs = [os.path.expanduser("~/.local/lib")]
library_dirs = [os.path.expanduser("~breddels/stow/64/lib")]
libraries = ["boost_python", "gsl", "gslcblas", "nlopt"]
defines = []
#if "darwin" in platform.system().lower():
extra_compile_args = ["-mfpmath=sse", "-O3", "-funroll-loops"]
#else:
#	extra_compile_args = ["-mfpmath=sse", "-msse4", "-Ofast", "-flto", "-march=native", "-funroll-loops"]
	#extra_compile_args = ["-mfpmath=sse", "-O3", "-funroll-loops"]
	#extra_compile_args = ["-mfpmath=sse", "-msse4a", "-O3", "-funroll-loops"]
extra_compile_args.extend(["-std=c++0x"])

include_dirs.append(os.path.join(get_python_inc(plat_specific=1), "numpy"))
include_dirs.append(os.path.join(numdir, "core", "include"))
import pyublas
include_dirs.append(os.path.join(os.path.dirname(pyublas.__file__), "include"))
#include_dirs.append("../../3rd/eigen")
#include_dirs.append("/Users/users/breddels/stow/64/stow/eigen2.0.12/include")
include_dirs.append("/Users/users/breddels/stow/64/include")

sources = "coordinate_systems.cpp df.cpp profile.cpp profile_2d.cpp sparse_grid_dynamic.cpp sparse_grid.cpp profile_axi.cpp profile3d.cpp transformation.cpp mesh2.cpp  poisson_fem2.cpp polynomial.cpp mesh.cpp functions.cpp galaxy.cpp aperture.cpp jeans.cpp gdfast.cpp orbit_integrator.cpp torus_spherical.cpp".split()

sources = [os.path.join("src", source) for source in sources]
print sources
extensions = [
	Extension("mab.gd.gdfast", sources,
                include_dirs=include_dirs,
                library_dirs=library_dirs,
                libraries=libraries,
                define_macros=defines,
                extra_compile_args=extra_compile_args
                )
]

setup(
    name="mab",
    version="0.2",
    data_files=DATA_FILES,
    ext_modules=extensions,
    package_data={'gavi': ['gavi/icons/*.png']},
)

